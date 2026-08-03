#!/usr/bin/env python3
"""Rebuild the deterministic real-data demo from the bundled EpiGym tables."""

from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
N_DOMAINS = 4
ROWS_PER_DOMAIN = 120


def main():
    reference = pd.read_csv(ROOT / "data" / "tsuboyama_stability_doubles.csv")
    predictions = pd.read_csv(
        ROOT / "baselines" / "tsuboyama_epistasis_predictions_per_double.csv",
        usecols=["DMS_id", "mutant", "eps_ProSST"],
    )

    merged = reference.merge(predictions, on=["DMS_id", "mutant"], how="inner")
    merged = merged.dropna(subset=["score", "additive_pred", "epistasis", "eps_ProSST"])

    counts = merged.groupby("DMS_id").size()
    domains = sorted(counts[counts >= ROWS_PER_DOMAIN].index)[:N_DOMAINS]
    if len(domains) != N_DOMAINS:
        raise RuntimeError("The bundled data do not contain four eligible demo domains")

    subset = pd.concat(
        [
            merged[merged["DMS_id"] == domain]
            .sort_values("mutant")
            .head(ROWS_PER_DOMAIN)
            for domain in domains
        ],
        ignore_index=True,
    )

    subset[["DMS_id", "mutant", "score", "additive_pred", "epistasis"]].to_csv(
        HERE / "demo_reference.csv", index=False, float_format="%.10g"
    )
    subset[["DMS_id", "mutant", "eps_ProSST"]].rename(
        columns={"eps_ProSST": "prediction"}
    ).to_csv(HERE / "demo_predictions.csv", index=False, float_format="%.10g")

    print(f"Wrote {len(subset)} real double mutants from {len(domains)} domains")
    for domain in domains:
        print(f"  {domain}: {ROWS_PER_DOMAIN}")


if __name__ == "__main__":
    main()
