# Is the EpiGym target real signal, or measurement noise?

The central target of EpiGym is the **epistatic residual**
`ε = effect(multi) − Σ effect(singles)`. Because ε is a *difference* of measurements, a
fair reviewer will ask: **is the ~0.32 prediction ceiling a genuine model limitation, or are
models simply failing to predict measurement noise?** This document collects the evidence
that ε carries real, learnable physical signal. Every numbered check is **reproducible from
the shipped files**:

```bash
python eval/validity_checks.py
```

**Falsification frame.** If ε were dominated by noise, then (1) it would be flat in 3-D
space, (2) no model would do better on contacting than distant pairs, (3) it would not sort
into physically-defined classes, (4) it would not agree across independent assays, and
(5) it would show no consistent scaling behaviour. **All five fail for the noise hypothesis.**

---

### 1 · ε is spatially structured (decays with residue distance)
Measured |ε| (per-domain z) falls monotonically with the Cβ–Cβ distance between the two
mutated residues — exactly what physical coupling predicts and what noise cannot produce:

| Cβ–Cβ distance | mean \|ε\| (z) | n |
|---|---|---|
| <8 Å (contact) | **1.46** | 38,450 |
| 8–10 Å | 1.36 | 7,894 |
| 10–13 Å | 0.76 | 3,018 |
| ≥13 Å | **0.45** | 1,450 |

76% of assayed doubles are in contact. Noise would be flat across these bins.

### 2 · Models predict contact ε better than distant ε
Per-domain |Spearman| on measured ε, split by distance (unsaturated subset, canonical T2):

| model | contact (<8 Å) | distant (≥8 Å) | contact better? |
|---|---|---|---|
| ProSST | **0.324** | 0.231 | yes |
| ESM2-650M | 0.257 | 0.224 | yes |
| EVmutation | 0.227 | 0.201 | yes |
| ESM1v | 0.190 | 0.152 | yes |
| ESM2-3B | 0.215 | 0.183 | yes |
| SPURS (structure-ΔΔG) | 0.085 | 0.097 | no* |

5 of 6 models capture *more* signal exactly where coupling is physically strongest.
\*SPURS is the lone exception precisely because it sits at the noise floor (|ρ|≈0.09)
everywhere — it captures essentially no ε signal at any distance, which is consistent with,
not counter to, the signal interpretation.

### 3 · ε sorts into physical classes with monotone magnitude (SKEMPI)
Classifying double-mutant ε from signed component ΔΔGs yields interpretable categories whose
magnitude increases exactly as the biophysics dictates — noise would not separate this way:

| type | share | mean \|ε\| (kcal/mol) |
|---|---|---|
| additive (negligible) | 33.7% | 0.15 |
| magnitude | 55.7% | 0.99 |
| sign | 8.7% | 1.93 |
| reciprocal-sign | 1.9% | 2.75 |

The ~11% sign/reciprocal-sign fraction is **unrepresentable by additive models in principle**
and carries the largest |ε|.

### 4 · The blind spot is consistent across four independent assay platforms
Stability (Tsuboyama), binding (SKEMPI), fitness (ProteinGym) and antibody escape (MAAD) use
unrelated measurement technologies with unrelated noise. The failure is the *same direction*
in all of them — e.g. on MAAD, modelling interactions adds **+0.10** AUROC when test variants
recombine *seen* mutations but **−0.01** when they introduce *novel* ones. Independent noise
sources do not align.

### 5 · ε prediction has a reproducible scaling signature *(paper; ED Fig. 2)*
Across the ESM2 family (8M→3B), epistasis-prediction accuracy traces a reproducible
inverted-U — it peaks near 150M parameters and then *declines* — while single-mutant accuracy
plateaus. A reproducible, non-monotone capacity dependence is a property of structured signal,
not of noise. *(This check uses per-model-size predictions that are not bundled in EpiGym v1;
see the paper.)*

---

### Protocol safeguards already built in
- **Unsaturated subset** (canonical T2): per domain, scoring is restricted to variants whose
  `score` and `additive_pred` lie within the 15–85th percentiles, removing the assay
  ceiling/floor where ΔG clipping would *manufacture* spurious ε.
- **Per-domain z-normalisation** of |ε| (`epistasis_z`) and a **≥15 doubles/domain** minimum
  guard against scale and small-sample artifacts.

### Rebuttal-ready summary
> The epistatic residual is not measurement noise. It is spatially organised (|ε| decays from
> 1.46 to 0.45 z with Cβ–Cβ distance; 76% of doubles are contacts), it is differentially
> *predictable* exactly where coupling is strongest (contact > distant for 5/6 models), it
> decomposes into physically-defined sign/magnitude classes with monotonically increasing
> magnitude (0.15→2.75 kcal/mol), it fails in the same direction across four independent assay
> platforms, and it has a reproducible inverted-U scaling signature. The ~0.32 ceiling
> therefore reflects a genuine, shared model limitation on real biophysical signal — which is
> exactly what makes it a worthwhile benchmark target.
