# Changelog

## v1.2 package build (2026-08-03)

Added software documentation and demonstration assets without changing benchmark values,
frozen splits, leaderboard values, or reference predictions:

- `requirements.txt` and the verified `requirements-tested.txt` environment;
- a deterministic four-domain Demo using real Tsuboyama measurements and ProSST predictions under `demo/`;
- an optional `--reference` argument so the evaluator can run on a small or external
  reference table through the same scoring path;
- strict prediction-key, coverage, numeric, and finite-value validation to prevent
  accidental partial or duplicated evaluations;
- an explicit `--allow-partial` mode that reports coverage for models with unavoidable
  scoring gaps, while complete coverage remains the default;
- a `--splits` option for evaluating a selected subset of MAAD split columns;
- updated the shipped construct-validity check to include all eight T2 reference models;
- explicit system, installation, runtime, usage, and reproduction instructions in
  `README.md`.
- excluded the upstream MAAD source table from the public archive; T5 now documents the
  official download page and requires a user-supplied reference table via `--reference`.

## v1.2 (2026-06-08)
Renamed the benchmark from **EpiBench** to **EpiGym** to avoid a name collision with unrelated tools already named "EpiBench". No changes to data, splits, leaderboards, the evaluator, or any results; internal version strings (left at 1.0 through v1.1) are now synchronized to 1.2.

## v1.1 (2026-06-05)
Added ESM-C, VenusREM, and PoET to the reference leaderboards.

- Generalization leaderboard (T1): 29 to 31 models; added zero-shot ESM-C
  (`esmc_600m`) and VenusREM.
- Epistasis-prediction leaderboard (T2): 6 to 8 models; added ESM-C and PoET.
  Existing model values were unchanged (ProSST 0.3221, ESM2-650M 0.3195,
  SPURS 0.0949; n=50).
- `tsuboyama_epistasis_predictions_per_double.csv`: added `eps_ESM-C` and `eps_PoET` columns.
- The fitted eight-model combined value changed from 0.46 to 0.52. This fitted value is
  not a leaderboard entry.
- README / DATASHEET / MANIFEST counts updated; `CHECKSUMS.sha256` regenerated.

## v1.0 (2026-06-04)
Initial EpiGym release.

- Four assay types with directly-measured epistasis: ProteinGym DMS (1.50M variants),
  Tsuboyama stability (50,812 doubles / 50 domains), SKEMPI binding (685 multi-mutants;
  621 typed doubles), MAAD antibody escape (46,816 rows).
- Frozen deployment splits (7 for ProteinGym, 10 for MAAD).
- Auxiliary annotations: Cβ–Cβ contact distance (Tsuboyama), epistasis-type labels (SKEMPI).
- Reference leaderboards from 29 models / 6 paradigms (generalization, epistasis prediction,
  antibody escape) + per-double reference predictions.
- `eval/evaluate.py` reference scorer; leaderboards regenerate from the same code path.
- `VALIDITY.md` + `eval/validity_checks.py`: descriptive checks for spatial distance,
  contact-stratified prediction, SKEMPI classes, and cross-dataset reference values.
- Machine-readable task cards (`tasks.json`), `MANIFEST.json`, `CHECKSUMS.sha256`.

### Distribution note
The MAAD source table is excluded from the archive. Its official download location is
listed in `README.md`.
