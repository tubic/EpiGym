# Changelog

## v1.2 (2026-06-08)
Renamed the benchmark from **EpiBench** to **EpiGym** to avoid a name collision with unrelated tools already named "EpiBench". No changes to data, splits, leaderboards, the evaluator, or any results; internal version strings (left at 1.0 through v1.1) are now synchronized to 1.2.

## v1.1 (2026-06-05)
Added three current state-of-the-art models to the reference leaderboards; the cross-paradigm epistasis blind spot is unchanged.

- Generalization leaderboard (T1): 29 → **31 models** — added zero-shot **ESM-C** (esmc_600m) and **VenusREM** (current ProteinGym leader); both stay flat across splits (no collapse).
- Epistasis-prediction leaderboard (T2): 6 → **8 models** — added **ESM-C** and the autoregressive, MSA-conditioned **PoET**; both median |ρ| ≈ 0.30–0.32, squarely in the blind spot. Existing models reproduce exactly (ProSST 0.3221, ESM2-650M 0.3195, SPURS 0.0949; n=50).
- `tsuboyama_epistasis_predictions_per_double.csv`: added `eps_ESM-C` and `eps_PoET` columns.
- Fitted all-model ceiling updated 0.46 → **0.52** (8 models; ~73% of ε variance still unexplained) — oracle context only, not a leaderboard entry.
- README / DATASHEET / MANIFEST counts updated; `CHECKSUMS.sha256` regenerated.

## v1.0 (2026-06-04)
Initial release accompanying *"Protein variant-effect models predict single-mutation effects but not epistasis."*

- Four assay types with directly-measured epistasis: ProteinGym DMS (1.50M variants),
  Tsuboyama stability (50,812 doubles / 50 domains), SKEMPI binding (685 multi-mutants;
  621 typed doubles), MAAD antibody escape (46,816 rows).
- Frozen deployment splits (7 for ProteinGym, 9 for MAAD).
- Auxiliary annotations: Cβ–Cβ contact distance (Tsuboyama), epistasis-type labels (SKEMPI).
- Reference leaderboards from 29 models / 6 paradigms (generalization, epistasis prediction,
  antibody escape) + per-double reference predictions.
- `eval/evaluate.py` reference scorer; leaderboards regenerate from the same code path.
- `VALIDITY.md` + `eval/validity_checks.py`: construct-validity evidence that ε is real
  signal, not measurement noise (spatial gradient, contact predictability, physical
  sign/magnitude classes, cross-assay consistency), reproducible from the shipped files.
- Machine-readable task cards (`tasks.json`), `MANIFEST.json`, `CHECKSUMS.sha256`.

### Known follow-ups before public release
- Confirm redistribution rights / final license per upstream dataset.
- Mint a DOI (e.g. Zenodo) and freeze.
