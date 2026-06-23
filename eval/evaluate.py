#!/usr/bin/env python3
"""Reference evaluator for EpiGym (v1.2).

Dependency-light: pandas, numpy, scipy only.

Score your own model
--------------------
  python eval/evaluate.py --task T2_epistasis_prediction --pred my_preds.csv
  python eval/evaluate.py --task T1_generalization --pred my_preds.csv --split split_distance_controlled
  python eval/evaluate.py --task T5_cross_assay_escape --pred my_preds.csv

Prediction CSV schema (one prediction per evaluated row):
  T1_generalization        : columns [assay_id, mutations, prediction]   (prediction = predicted fitness)
  T2_epistasis_prediction  : columns [DMS_id, mutant, prediction]        (prediction = predicted epistasis residual)
  T5_cross_assay_escape    : columns [row_id, prediction]                (prediction = escape score / probability)

Reproduce the reference numbers
-------------------------------
  python eval/evaluate.py --selftest
"""
import argparse, json, os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def _load(rel):
    return pd.read_csv(os.path.join(ROOT, rel))


def auroc(y, s):
    """Rank-based AUROC (no sklearn dependency)."""
    y = np.asarray(y); s = np.asarray(s)
    pos, neg = (y == 1), (y == 0)
    npos, nneg = pos.sum(), neg.sum()
    if npos == 0 or nneg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s)); ranks[order] = np.arange(1, len(s) + 1)
    # average ranks for ties
    _, inv, cnt = np.unique(s, return_inverse=True, return_counts=True)
    csum = np.cumsum(cnt); start = csum - cnt
    avg = (start + csum + 1) / 2.0
    ranks = avg[inv]
    return float((ranks[pos].sum() - npos * (npos + 1) / 2) / (npos * nneg))


def eval_T1(pred, split="split_distance_controlled"):
    d = _load("data/proteingym_dms.csv.gz")
    m = d.merge(pred, on=["assay_id", "mutations"], how="inner")
    test = m[m[split] == "test"]
    rs = []
    for a, g in test.groupby("assay_id"):
        g = g.dropna(subset=["score", "prediction"])
        if len(g) >= 10:
            r = spearmanr(g.score, g.prediction).correlation
            if r == r:
                rs.append(r)
    return dict(task="T1_generalization", split=split,
                mean_spearman=round(float(np.mean(rs)), 4), n_assays=len(rs))


def eval_T2(pred, unsaturated=True):
    d = _load("data/tsuboyama_stability_doubles.csv")
    m = d.merge(pred, on=["DMS_id", "mutant"], how="inner")
    rs = []
    for dom, g in m.groupby("DMS_id"):
        if unsaturated:
            # canonical protocol: restrict to the unsaturated core so the stability
            # measurement ceiling/floor does not manufacture spurious epistasis
            slo, shi = g.score.quantile([.15, .85])
            alo, ahi = g.additive_pred.quantile([.15, .85])
            g = g[(g.score > slo) & (g.score < shi) & (g.additive_pred > alo) & (g.additive_pred < ahi)]
        g = g.dropna(subset=["epistasis", "prediction"])
        if len(g) >= 15:
            r = spearmanr(g.epistasis, g.prediction).correlation
            if r == r:
                rs.append(abs(r))
    return dict(task="T2_epistasis_prediction",
                protocol="unsaturated subset [15,85]%" if unsaturated else "all doubles",
                median_abs_spearman=round(float(np.median(rs)), 4),
                mean_abs_spearman=round(float(np.mean(rs)), 4), n_domains=len(rs))


def eval_T5(pred, splits=None):
    d = _load("data/maad_escape.csv.gz")
    m = d.merge(pred, on="row_id", how="inner")
    cols = splits or [c for c in d.columns if c.startswith("split_")]
    out = {}
    for s in cols:
        t = m[m[s] == "test"].dropna(subset=["label", "prediction"])
        if t["label"].nunique() == 2:
            out[s] = round(auroc(t["label"].values, t["prediction"].values), 4)
    return dict(task="T5_cross_assay_escape", auroc_by_split=out)


def selftest():
    """Reproduce reference numbers from the bundled baseline predictions."""
    print("Self-test: scoring bundled reference predictions")
    perd = _load("baselines/tsuboyama_epistasis_predictions_per_double.csv")
    for model in ["ProSST", "ESM2-650M", "SPURS"]:
        pred = perd[["DMS_id", "mutant"]].copy()
        pred["prediction"] = perd[f"eps_{model}"]
        r = eval_T2(pred)
        print(f"  T2 {model:10s}: median|rho|={r['median_abs_spearman']:.4f}  (n={r['n_domains']})")
    print("Expected (baselines/tsuboyama_epistasis_leaderboard.csv): ProSST~0.32, ESM2-650M~0.32, SPURS~0.09")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=["T1_generalization", "T2_epistasis_prediction", "T5_cross_assay_escape"])
    ap.add_argument("--pred", help="predictions CSV")
    ap.add_argument("--split", default="split_distance_controlled")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    else:
        pred = pd.read_csv(a.pred)
        fn = {"T1_generalization": lambda p: eval_T1(p, a.split),
              "T2_epistasis_prediction": eval_T2,
              "T5_cross_assay_escape": eval_T5}[a.task]
        print(json.dumps(fn(pred), indent=2))
