# Datasheet — EpiGym v1.2

Following *Datasheets for Datasets* (Gebru et al.).

## Purpose
EpiGym measures compositional generalization and prediction of the non-additive component
of multi-mutation effects.

## Composition
The package covers four assay types. ProteinGym, Tsuboyama, and SKEMPI tables include the
epistatic residual `ε = effect(multi) − Σ effect(singles)`; MAAD uses a binary escape
label.

| File | Unit / phenotype | Rows | ε target | Splits |
|---|---|---|---|---|
| `data/proteingym_dms.csv.gz` | DMS fitness (assay-specific) | 1,504,667 variants (≈99% multi-mutant) | `epistasis`, `abs_epistasis` | 7 (`split_*`) |
| `data/proteingym_assays.csv` | wildtype sequences + metadata | 12 assays | — | — |
| `data/tsuboyama_stability_doubles.csv` | folding ΔG (kcal/mol) | 50,812 doubles, 50 domains | `epistasis` (+`epistasis_z`) | — (zero-shot task) |
| `data/tsuboyama_all_variants.csv.gz` | folding ΔG + mutant sequence | 117,811 singles+doubles | — | — |
| `data/skempi_binding.csv` | binding ΔΔG (kcal/mol) | 685 multi-mutants (k=2–8) | `epistasis` | — |
| `data/skempi_epistasis_type.csv` | signed component ΔΔGs | 621 doubles | `epistasis` + `type` | — |
| user-supplied MAAD reference table (not bundled) | antibody×variant neutralization | 46,816 (label 0/1) | `label` (escape) | 10 (`split_*`) |

`epistasis_type` ∈ {additive, magnitude, sign, reciprocal-sign} (≈34/56/9/2%). The
sign and reciprocal-sign classes account for about 11% of the typed rows.

## Collection / provenance
EpiGym reuses and annotates the following public datasets:

- **ProteinGym** — Notin et al., *NeurIPS* 2023. Source DMS assays compiled by ProteinGym.
- **Tsuboyama mega-scale stability** — Tsuboyama et al., *Nature* 2023.
- **SKEMPI 2.0** — Jankauskaite et al., *Bioinformatics* 2019.
- **MAAD antibody–variant neutralization** — Li et al., *Protein & Cell* 2026;
  obtain the source data from https://www.raabmd.org/raab/download.
- **AlphaFold structures** (for Cβ–Cβ distance) — used to annotate Tsuboyama domains.

EpiGym supplies the derived ε values, frozen evaluation splits, epistasis-type labels,
contact distances, and reference leaderboards.

## Preprocessing
- **ε / additive_pred**: for DMS and stability, `additive_pred = Σ single-mutation effects`
  measured in the *same* assay; `ε = measured(multi) − additive_pred`. `epistasis_z` is the
  per-domain z-normalised |ε| (Tsuboyama).
- **Splits** (per-row `train/test/val/unused`): `random`; `single_to_multi` (train singles →
  test multi); `low_to_high_k2` (train ≤k → test higher order); `seen_marginal_unseen_combo`
  (test novel *combinations* of seen single positions); `epistasis_enriched`; and
  `distance_controlled` (test spatially-coupled pairs). MAAD: `random`, antibody-group /
  variant holdouts, and viral temporal splits (`pre_omicron_to_omicron`, `ba12_to_xbb_bq`),
  `low_to_high_burden`, `cross_virus_hardest`.
- **Cβ–Cβ distance**: between the two mutated residues from AlphaFold structures (Gly→Cα).
- **SKEMPI type**: from signed component ΔΔGs with a 0.3 kcal/mol threshold.

## Recommended uses
T1–T5 are defined in `tasks.json`. The T2 protocol restricts each Tsuboyama domain to rows
whose `score` and `additive_pred` are within the 15–85th percentiles before scoring. This
filter reduces the effect of assay saturation. `eval/evaluate.py` and the bundled
leaderboard use the same filter.

Do not compare raw `score` or `ε` values across assays because their units and scales
differ. Do not treat `epistasis_type` counts and `skempi_binding` rows as the same sample;
they are independent extractions with different row counts.

## Distribution & license
Derived annotations (ε, splits, type labels, distances, leaderboards) — **CC-BY-4.0**.
Underlying measurements/sequences remain under their original sources' terms. **Before any
public/DOI redistribution, confirm you may redistribute each upstream dataset**; some may
require linking to the source rather than re-hosting raw data.

The MAAD source table is therefore not included in this archive. Users must obtain it from
the official MAAD download page and provide an EpiGym-compatible T5 reference table using
the evaluator's `--reference` option. The derived MAAD leaderboard is included.

## Maintenance
Versioned (`v1.2`) and integrity-checked via `CHECKSUMS.sha256`. The archive is a frozen
evaluation release; the upstream-dependent assembly script is maintained in the project
repository rather than bundled here. Report issues or submit models via
https://github.com/tubic/EpiGym.

## Known limitations
- DMS `score` is assay-relative; cross-assay aggregation must use rank metrics per assay.
- Stability ΔG saturates at the assay limits → T2 uses the unsaturated subset.
- SKEMPI multi-mutant coverage is modest (685 rows; 621 typed doubles, 31 repeated (pdb,muts)).
- MAAD labels are binarised neutralization with weak-label provenance.
