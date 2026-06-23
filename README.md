# EpiGym — v1.2

*A diagnostic benchmark for compositional (multi-mutation) generalization and the direct
prediction of measured epistasis in protein variant-effect models.*

This release packages everything needed to **reproduce its evaluation and score new models** on the same footing. EpiGym is a **data + reference-leaderboard + evaluator** release; it does **not** ship model weights or embeddings.

---

## Why EpiGym

Standard variant-effect leaderboards reward predicting *single*-mutation effects and
average over deployment regimes, so they hide a specific failure: models do not predict
the **epistatic residual** — the part of a multi-mutant's effect that is *not* the sum of
its single mutations, `ε = effect(AB) − [effect(A) + effect(B)]`. EpiGym isolates that
failure with (i) **deployment-like generalization splits**, (ii) **directly-measured ε** as
an explicit prediction target, and (iii) coverage across **four assay types** so the result
is not an artifact of one phenotype.

**Headline reference results (this release):**

| Task | Best reference | Additive floor | Note |
|---|---|---|---|
| Predict measured stability epistasis (T2) | ProSST **|ρ|≈0.32** | ≈0 by construction | every paradigm ≤0.32 (incl. the latest ESM-C / PoET); a *fitted* 8-model combination reaches only ≈**0.52** → **~73% of ε variance is unexplained** |
| Structure-ΔΔG on the same task | SPURS **|ρ|≈0.09** | — | the most structure-aware model is the *worst* at epistasis |
| High-order generalization (T1) | see leaderboard | — | every paradigm degrades from `random` to `distance_controlled`; coevolution degrades least |
| Antibody escape (T5) | interaction model | additive | interactions help only when test variants *recombine seen* mutations (+0.09–0.11 AUROC), not when they introduce *novel* ones (≈0 / negative) |

The **0.52 fitted-combination ceiling** is an oracle (it uses the target to fit) and is
reported here as *context only* — it is **not** a leaderboard entry.

> **Is the target real signal or noise?** Because ε is a difference of measurements, EpiGym
> ships an explicit construct-validity argument with reproducible checks — see `VALIDITY.md`
> (`python eval/validity_checks.py`): ε is spatially structured (decays with residue
> distance), differentially predictable at contacts, sorts into physical sign/magnitude
> classes, and is consistent across all four assays.

---

## Contents

```
EpiGym_v1/
├── README.md                 ·  this file
├── DATASHEET.md              ·  provenance, licenses, composition, caveats
├── VALIDITY.md               ·  evidence that ε is real signal, not noise (+ reproducible checks)
├── CHANGELOG.md
├── tasks.json                ·  machine-readable task cards
├── MANIFEST.json             ·  per-file row/col/size/notes
├── CHECKSUMS.sha256
├── data/
│   ├── proteingym_dms.csv.gz            ·  1.50M variants · fitness · additive_pred · ε · 7 splits
│   ├── proteingym_assays.csv            ·  12 assays · wildtype sequence + metadata
│   ├── tsuboyama_stability_doubles.csv  ·  50,812 doubles · ddG · ε (raw/|.|/z) · Cβ–Cβ distance
│   ├── tsuboyama_all_variants.csv.gz    ·  117,811 singles+doubles · ddG + mutant sequence
│   ├── skempi_binding.csv               ·  685 multi-mutants · binding ΔΔG · ε · order k
│   ├── skempi_epistasis_type.csv        ·  621 doubles · signed component ΔΔGs · epistasis type
│   └── maad_escape.csv.gz               ·  46,816 antibody×variant rows · label · 9 splits
├── baselines/
│   ├── proteingym_generalization_leaderboard.csv      ·  31 models × 12 assays × 7 splits
│   ├── tsuboyama_epistasis_leaderboard.csv            ·  8 models · median |ρ| on measured ε
│   ├── tsuboyama_epistasis_predictions_per_double.csv ·  per-double measured + predicted ε
│   └── maad_escape_leaderboard.csv                    ·  additive vs interaction AUROC per split
└── eval/
    ├── evaluate.py           ·  reference scorer (pandas + numpy + scipy)
    └── validity_checks.py    ·  construct-validity checks (ε is signal, not noise)
```

Reconstruct any ProteinGym mutant sequence from `proteingym_assays.wildtype_sequence` +
the `mutations` string (e.g. `A123C:D200E`, 1-indexed).

---

## Tasks

Full definitions in `tasks.json`. Summary:

| ID | Task | Data | Target | Metric |
|---|---|---|---|---|
| **T1** | High-order generalization | `proteingym_dms` | `score` | Spearman ρ per assay, mean; report drop `random`→`distance_controlled` |
| **T2** | **Direct epistasis prediction** (blind-spot task) | `tsuboyama_stability_doubles` | `epistasis` | per-domain \|Spearman\| on the unsaturated subset, median over domains |
| **T3** | Epistasis-stratified accuracy | `proteingym_dms` | `score` | within-quintile Spearman by `abs_epistasis`; the top-quintile "cliff" |
| **T4** | Contact-resolved epistasis | `tsuboyama_stability_doubles` | `epistasis` | \|Spearman\| at contacts (<8 Å) vs distant (≥8 Å) |
| **T5** | Antibody-escape generalization | `maad_escape` | `label` | AUROC per split |
| **T*** | Binding epistasis (auxiliary) | `skempi_binding` | `epistasis` | Spearman by order *k*; type composition |

**Splits** are frozen, stored as per-row columns (`split_*`) with values
`train` / `test` / `val` / `unused`. ProteinGym splits: `random`, `single_to_multi`,
`low_to_high_k1` (≡ `single_to_multi` by construction), `low_to_high_k2`,
`seen_marginal_unseen_combo`, `epistasis_enriched`, `distance_controlled` (defined for the
6 assays with structures). MAAD adds viral temporal / antibody-holdout splits (see
`tasks.json`).

---

## How to evaluate

```bash
# reproduce the reference numbers end-to-end
python eval/evaluate.py --selftest
#   T2 ProSST : median|rho|=0.3221 (n=50)   ← matches tsuboyama_epistasis_leaderboard.csv

# reproduce the construct-validity evidence (ε is signal, not noise)
python eval/validity_checks.py

# score your own model
python eval/evaluate.py --task T2_epistasis_prediction --pred my_preds.csv
python eval/evaluate.py --task T1_generalization       --pred my_preds.csv --split split_distance_controlled
python eval/evaluate.py --task T5_cross_assay_escape   --pred my_preds.csv
```

Prediction CSV schema (one row per evaluated variant):

| Task | Required columns |
|---|---|
| T1 | `assay_id`, `mutations`, `prediction` (predicted fitness) |
| T2 | `DMS_id`, `mutant`, `prediction` (predicted **epistasis residual**) |
| T5 | `row_id`, `prediction` (escape score / probability) |

The evaluator and the shipped leaderboards use **identical** code paths, so your numbers are
directly comparable to the reference baselines.

---

## Citation

```bibtex
@article{epistasis_blindspot_2026,
  title  = {Protein variant-effect models predict single-mutation effects but not epistasis},
  author = {<authors>},
  year   = {2026},
  note   = {EpiGym v1.2}
}
```

Please **also cite the upstream data sources** (ProteinGym, Tsuboyama et al. 2023,
SKEMPI 2.0, and the MAAD compilation) listed in `DATASHEET.md`.

## License & caveats

Derived annotations in EpiGym (measured ε, frozen splits, epistasis-type labels,
Cβ–Cβ distances, leaderboards) are released under **CC-BY-4.0**. The underlying assay
measurements and sequences remain subject to their original sources' licenses — see
`DATASHEET.md`. **Verify redistribution rights for each upstream dataset before any public
DOI release.**
