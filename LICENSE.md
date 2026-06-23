# License & Notice

## Derived annotations — CC-BY-4.0
The contributions original to this benchmark are licensed **Creative Commons
Attribution 4.0 International (CC-BY-4.0)**:

- the directly-measured epistasis values (`epistasis`, `abs_epistasis`, `epistasis_z`)
  and additive predictions (`additive_pred`),
- the frozen generalization splits (`split_*` columns),
- the SKEMPI epistasis-type labels, the Cβ–Cβ contact distances,
- the reference leaderboards, task cards, and evaluation code.

Attribution: cite *"Protein variant-effect models predict single-mutation effects but not epistasis"* (2026) and EpiGym v1.2.

## Underlying data — original terms apply
The benchmark **re-uses** measurements and sequences from third-party datasets that retain
their own licenses and citation requirements:

| Source | Reference | Terms |
|---|---|---|
| ProteinGym (DMS assays) | Notin et al. 2023 | per ProteinGym + each original DMS study |
| Tsuboyama mega-scale stability | Tsuboyama et al. 2023 | per the original release |
| SKEMPI 2.0 | Jankauskaite et al. 2019 | per SKEMPI |
| MAAD antibody–variant neutralization | community compilation | per source |
| AlphaFold structures (distance annotation) | Jumper et al. 2021 | CC-BY-4.0 |

**Before any public redistribution or DOI minting, verify that you may re-host the raw
upstream measurements/sequences.** If a source does not permit re-hosting, distribute the
derived annotations + split indices and a script that joins them to the user-downloaded
source data.
