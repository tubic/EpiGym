#!/usr/bin/env python3
"""Run and verify the small real-data EpiGym T2 demonstration."""

import json
import sys
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "eval"))

from evaluate import eval_T2  # noqa: E402


def main():
    predictions = pd.read_csv(HERE / "demo_predictions.csv")
    result = eval_T2(predictions, reference=HERE / "demo_reference.csv")

    with (HERE / "expected_output.json").open(encoding="ascii") as handle:
        expected = json.load(handle)

    print(json.dumps(result, indent=2))
    if result != expected:
        raise SystemExit(f"DEMO FAILED: expected {expected}, received {result}")
    print("DEMO PASS: output matches demo/expected_output.json")


if __name__ == "__main__":
    main()
