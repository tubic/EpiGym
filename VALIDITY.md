# Validation checks

The EpiGym T2 target is the epistatic residual
`ε = effect(multi) − Σ effect(singles)`. This document records four descriptive checks
implemented in `eval/validity_checks.py`.

Run the checks from the package root:

```bash
python eval/validity_checks.py
```

## 1. Residue-distance stratification

Measured absolute epistasis, normalized within each domain, is summarized by the Cβ–Cβ
distance between the two mutated residues.

| Cβ–Cβ distance | mean \|ε\| (z) | n |
|---|---:|---:|
| <8 Å | 1.46 | 38,450 |
| 8–10 Å | 1.36 | 7,894 |
| 10–13 Å | 0.76 | 3,018 |
| ≥13 Å | 0.45 | 1,450 |

The mean decreases across the four distance bins. Contact pairs account for 76% of the
double-mutant rows in this table.

## 2. Contact and distant subsets

The table reports per-domain absolute Spearman correlation on the unsaturated T2 subset,
calculated separately for contact and distant residue pairs.

| model | contact (<8 Å) | distant (≥8 Å) |
|---|---:|---:|
| ProSST | 0.324 | 0.231 |
| ESM2-650M | 0.257 | 0.224 |
| EVmutation | 0.227 | 0.201 |
| ESM1v | 0.190 | 0.152 |
| ESM2-3B | 0.215 | 0.183 |
| SPURS | 0.085 | 0.097 |
| ESM-C | 0.319 | 0.210 |
| PoET | 0.303 | 0.170 |

Seven of the eight reference models have a larger value on the contact subset. SPURS has
0.085 on the contact subset and 0.097 on the distant subset.

## 3. SKEMPI epistasis classes

Double mutants are classified from the signed component ΔΔG values. The table reports the
fraction of typed rows and the mean absolute epistatic residual.

| type | share | mean \|ε\| (kcal/mol) |
|---|---:|---:|
| additive | 33.7% | 0.15 |
| magnitude | 55.7% | 0.99 |
| sign | 8.7% | 1.93 |
| reciprocal-sign | 1.9% | 2.75 |

## 4. Cross-dataset reference summary

The script reads the bundled Tsuboyama and MAAD leaderboards. It reports the highest
Tsuboyama T2 correlation and the mean MAAD interaction-model AUROC gain for the configured
recombination and novel-mutation split groups. This check uses the derived MAAD leaderboard;
the MAAD source table is not required.

## Evaluation controls

- T2 scoring uses rows whose `score` and `additive_pred` are within the 15th and 85th
  percentiles for each domain.
- Absolute epistasis is normalized within each domain for the distance summary.
- Per-domain correlations require at least 15 double mutants.
