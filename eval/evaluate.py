#!/usr/bin/env python3
"""Reference evaluator for EpiGym (v1.2).

Dependency-light: pandas, numpy, scipy only.

Score your own model
--------------------
  python eval/evaluate.py --task T2_epistasis_prediction --pred my_preds.csv
  python eval/evaluate.py --task T1_generalization --pred my_preds.csv --split split_distance_controlled
  python eval/evaluate.py --task T5_cross_assay_escape \
      --reference path/to/maad_escape.csv.gz --pred my_preds.csv
  python eval/evaluate.py --task T5_cross_assay_escape --splits split_random \
      --reference path/to/maad_escape.csv.gz --pred my_preds.csv

Run against a small or external reference table
------------------------------------------------
  python eval/evaluate.py --task T2_epistasis_prediction \
      --reference demo/demo_reference.csv --pred demo/demo_predictions.csv

Prediction CSV schema (one prediction per evaluated row):
  T1_generalization        : columns [assay_id, mutations, prediction]   (prediction = predicted fitness)
  T2_epistasis_prediction  : columns [DMS_id, mutant, prediction]        (prediction = predicted epistasis residual)
  T5_cross_assay_escape    : columns [row_id, prediction]                (prediction = escape score / probability)

Reproduce the reference numbers
-------------------------------
  python eval/evaluate.py --selftest
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VERSION = "1.2"
T5_DEFAULT_SPLITS = [
    "split_random",
    "split_seen_positions_unseen_combinations",
    "split_ba12_to_xbb_bq",
    "split_low_to_high_burden",
    "split_pre_omicron_to_omicron",
    "split_cross_virus_hardest",
]


def _load(rel):
    return pd.read_csv(os.path.join(ROOT, rel))


def auroc(y, s):
    """Rank-based AUROC (no sklearn dependency)."""
    y = np.asarray(y)
    s = np.asarray(s)
    pos, neg = (y == 1), (y == 0)
    npos, nneg = pos.sum(), neg.sum()
    if npos == 0 or nneg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s))
    ranks[order] = np.arange(1, len(s) + 1)
    # average ranks for ties
    _, inv, cnt = np.unique(s, return_inverse=True, return_counts=True)
    csum = np.cumsum(cnt)
    start = csum - cnt
    avg = (start + csum + 1) / 2.0
    ranks = avg[inv]
    return float((ranks[pos].sum() - npos * (npos + 1) / 2) / (npos * nneg))


def _reference(default_rel, reference=None):
    """Load the bundled task data or an explicitly supplied reference CSV."""
    return pd.read_csv(reference) if reference else _load(default_rel)


def _require_columns(frame, columns, label):
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def _validate_predictions(
    pred, reference, keys, expected_mask=None, allow_partial=False
):
    """Validate keys and return prediction coverage for the evaluation subset."""
    _require_columns(pred, [*keys, "prediction"], "prediction file")

    reference_duplicates = int(reference.duplicated(keys).sum())
    if reference_duplicates:
        raise ValueError(
            f"reference file contains {reference_duplicates} duplicate task keys"
        )

    prediction_duplicates = int(pred.duplicated(keys).sum())
    if prediction_duplicates:
        raise ValueError(
            f"prediction file contains {prediction_duplicates} duplicate task keys"
        )

    expected = (
        reference.loc[expected_mask, keys]
        if expected_mask is not None
        else reference[keys]
    )
    if expected.empty:
        raise ValueError(
            "reference file contains no rows in the requested evaluation subset"
        )

    clean = pred.copy()
    raw_prediction = clean["prediction"]
    clean["prediction"] = pd.to_numeric(raw_prediction, errors="coerce")
    nonnumeric = raw_prediction.notna() & clean["prediction"].isna()
    if nonnumeric.any():
        raise ValueError(
            f"prediction file contains {int(nonnumeric.sum())} non-numeric predictions"
        )
    nonfinite = clean["prediction"].notna() & ~np.isfinite(clean["prediction"])
    if nonfinite.any():
        raise ValueError(
            f"prediction file contains {int(nonfinite.sum())} non-finite predictions"
        )

    expected_index = expected.set_index(keys).index
    supplied = clean.set_index(keys)["prediction"].reindex(expected_index)

    missing = int(supplied.isna().sum())
    if missing and not allow_partial:
        raise ValueError(
            "prediction coverage is incomplete for the requested evaluation subset: "
            f"{missing} predictions are missing; use --allow-partial only when the "
            "model cannot score those rows"
        )
    predicted = len(expected) - missing
    if predicted == 0:
        raise ValueError(
            "prediction file contains no usable rows in the evaluation subset"
        )

    coverage = {
        "n_expected": int(len(expected)),
        "n_predicted": int(predicted),
        "prediction_coverage": round(float(predicted / len(expected)), 4),
    }
    return clean, coverage


def _unsaturated_mask(frame):
    """Return the canonical per-domain T2 mask before predictions are joined."""
    keep = pd.Series(False, index=frame.index)
    for _, group in frame.groupby("DMS_id"):
        score_low, score_high = group.score.quantile([0.15, 0.85])
        additive_low, additive_high = group.additive_pred.quantile([0.15, 0.85])
        selected = (
            group.score.gt(score_low)
            & group.score.lt(score_high)
            & group.additive_pred.gt(additive_low)
            & group.additive_pred.lt(additive_high)
        )
        keep.loc[group.index] = selected
    return keep


def eval_T1(
    pred, split="split_distance_controlled", reference=None, allow_partial=False
):
    d = _reference("data/proteingym_dms.csv.gz", reference)
    _require_columns(d, ["assay_id", "mutations", "score", split], "reference file")
    keys = ["assay_id", "mutations"]
    pred, coverage = _validate_predictions(
        pred, d, keys, d[split].eq("test"), allow_partial
    )
    m = d.merge(pred, on=keys, how="inner", validate="one_to_one")
    test = m[m[split] == "test"]
    rs = []
    for a, g in test.groupby("assay_id"):
        g = g.dropna(subset=["score", "prediction"])
        if len(g) >= 10:
            r = spearmanr(g.score, g.prediction).correlation
            if r == r:
                rs.append(r)
    if not rs:
        raise ValueError(
            "no assay had at least 10 finite rows with a defined correlation"
        )
    return dict(
        task="T1_generalization",
        split=split,
        mean_spearman=round(float(np.mean(rs)), 4),
        n_assays=len(rs),
        **coverage,
    )


def eval_T2(pred, unsaturated=True, reference=None, allow_partial=False):
    d = _reference("data/tsuboyama_stability_doubles.csv", reference)
    _require_columns(
        d,
        ["DMS_id", "mutant", "score", "additive_pred", "epistasis"],
        "reference file",
    )
    keys = ["DMS_id", "mutant"]
    evaluation_mask = (
        _unsaturated_mask(d) if unsaturated else pd.Series(True, index=d.index)
    )
    pred, coverage = _validate_predictions(
        pred, d, keys, evaluation_mask, allow_partial
    )
    m = d.loc[evaluation_mask].merge(pred, on=keys, how="inner", validate="one_to_one")
    rs = []
    for dom, g in m.groupby("DMS_id"):
        g = g.dropna(subset=["epistasis", "prediction"])
        if len(g) >= 15:
            r = spearmanr(g.epistasis, g.prediction).correlation
            if r == r:
                rs.append(abs(r))
    if not rs:
        raise ValueError(
            "no domain had at least 15 finite rows with a defined correlation"
        )
    return dict(
        task="T2_epistasis_prediction",
        protocol="unsaturated subset [15,85]%" if unsaturated else "all doubles",
        median_abs_spearman=round(float(np.median(rs)), 4),
        mean_abs_spearman=round(float(np.mean(rs)), 4),
        n_domains=len(rs),
        **coverage,
    )


def eval_T5(pred, splits=None, reference=None, allow_partial=False):
    d = _reference("data/maad_escape.csv.gz", reference)
    _require_columns(d, ["row_id", "label"], "reference file")
    cols = splits or T5_DEFAULT_SPLITS
    if not cols:
        raise ValueError("reference file contains no split_* columns to evaluate")
    _require_columns(d, cols, "reference file")
    expected_mask = d[cols].eq("test").any(axis=1)
    pred, coverage = _validate_predictions(
        pred, d, ["row_id"], expected_mask, allow_partial
    )
    m = d.merge(pred, on="row_id", how="inner", validate="one_to_one")
    out = {}
    for s in cols:
        t = m[m[s] == "test"].dropna(subset=["label", "prediction"])
        if t["label"].nunique() == 2:
            out[s] = round(auroc(t["label"].values, t["prediction"].values), 4)
    if not out:
        raise ValueError("no requested split contains both outcome classes")
    return dict(task="T5_cross_assay_escape", auroc_by_split=out, **coverage)


def selftest():
    """Reproduce reference numbers from the bundled baseline predictions."""
    print("Self-test: scoring bundled reference predictions")
    perd = _load("baselines/tsuboyama_epistasis_predictions_per_double.csv")
    leaderboard = _load("baselines/tsuboyama_epistasis_leaderboard.csv")
    for row in leaderboard.itertuples(index=False):
        model = row.model
        pred = perd[["DMS_id", "mutant"]].copy()
        pred["prediction"] = perd[f"eps_{model}"]
        partial = bool(pred["prediction"].isna().any())
        result = eval_T2(pred, allow_partial=partial)
        if (
            result["median_abs_spearman"] != row.median_abs_rho
            or result["mean_abs_spearman"] != row.mean_abs_rho
            or result["n_domains"] != row.n_domains
        ):
            raise AssertionError(f"{model} does not match the shipped leaderboard")
        print(
            f"  T2 {model:10s}: median|rho|={result['median_abs_spearman']:.4f}  "
            f"(n={result['n_domains']}, coverage={result['prediction_coverage']:.1%})"
        )
    print("PASS: all eight models match baselines/tsuboyama_epistasis_leaderboard.csv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Score predictions with the EpiGym evaluator"
    )
    ap.add_argument(
        "--version", action="version", version=f"EpiGym evaluator {VERSION}"
    )
    ap.add_argument(
        "--task",
        choices=[
            "T1_generalization",
            "T2_epistasis_prediction",
            "T5_cross_assay_escape",
        ],
    )
    ap.add_argument("--pred", help="predictions CSV")
    ap.add_argument(
        "--split",
        default="split_distance_controlled",
        help="T1 split column (default: split_distance_controlled)",
    )
    ap.add_argument(
        "--splits",
        help="comma-separated T5 split columns; defaults to the six leaderboard splits",
    )
    ap.add_argument(
        "--reference",
        help=(
            "task reference CSV; optional for T1/T2 (bundled defaults), required for "
            "T5 because MAAD source data are not redistributed"
        ),
    )
    ap.add_argument(
        "--allow-partial",
        action="store_true",
        help="allow missing prediction rows and report coverage (default: require complete)",
    )
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    else:
        if not a.task:
            ap.error("--task is required unless --selftest is used")
        if not a.pred:
            ap.error("--pred is required unless --selftest is used")
        if a.splits and a.task != "T5_cross_assay_escape":
            ap.error("--splits can only be used with T5_cross_assay_escape")
        try:
            pred = pd.read_csv(a.pred)
        except (FileNotFoundError, OSError, pd.errors.ParserError) as exc:
            ap.error(f"could not read prediction file: {exc}")
        t5_splits = (
            [s.strip() for s in a.splits.split(",") if s.strip()] if a.splits else None
        )
        fn = {
            "T1_generalization": lambda p: eval_T1(
                p, a.split, a.reference, a.allow_partial
            ),
            "T2_epistasis_prediction": lambda p: eval_T2(
                p, reference=a.reference, allow_partial=a.allow_partial
            ),
            "T5_cross_assay_escape": lambda p: eval_T5(
                p,
                splits=t5_splits,
                reference=a.reference,
                allow_partial=a.allow_partial,
            ),
        }[a.task]
        try:
            result = fn(pred)
        except (FileNotFoundError, OSError, ValueError, KeyError) as exc:
            ap.error(str(exc))
        print(json.dumps(result, indent=2))
