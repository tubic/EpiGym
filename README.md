# EpiGym v1.2: software and demonstration package

EpiGym is an evaluation package for compositional (multi-mutation) generalization and
direct prediction of measured epistasis in protein variant-effect models. This package
contains the source code for the reference evaluator, the redistributable benchmark tables,
a small real-data demonstration dataset, and reference predictions. The upstream MAAD
source table is not redistributed in this archive; T5 users obtain it separately from the
official MAAD download page as described in section 4.

The evaluator is distributed as readable Python source code in `eval/` and `demo/`. It is
not a compiled binary and does not require a compiler, a GPU, or model weights. The bundled
Demo and T1/T2 tasks work offline after the Python dependencies have been installed; T5
requires a one-time download of the upstream MAAD data.

## Package guide

| Required item | Location in this package |
|---|---|
| Standalone source code | `eval/evaluate.py`, `eval/validity_checks.py`, `demo/*.py` |
| Small real demo dataset | `demo/demo_reference.csv`, `demo/demo_predictions.csv` |
| System requirements and tested versions | README section 1; `requirements.txt`; `requirements-tested.txt` |
| Installation instructions and install time | README section 2 |
| Demo command, output, and runtime | README section 3; `demo/expected_output.json` |
| Instructions for using your own data | README section 4 |
| Reproduction commands (optional item) | README section 5; `VALIDITY.md` |

## 1. System requirements

### Operating systems and Python

- CPython **3.10 or newer**.
- Linux, macOS, and Windows are suitable because the evaluator uses only portable Python,
  NumPy, pandas, and SciPy code. The package was tested on **Ubuntu 22.04.5 LTS
  (x86_64)**. macOS and Windows were not part of the verification run.
- A virtual environment is recommended. No system-wide installation is required.

### Runtime dependencies

The evaluator requires only the following runtime packages; compatible ranges are in
`requirements.txt`:

| Package | Supported version constraint |
|---|---|
| Python | 3.10+ |
| NumPy | >=1.24, <3 |
| pandas | >=2.0, <4 |
| SciPy | >=1.10, <2 |

pandas also installs its transitive runtime dependencies `python-dateutil>=2.8.2` and
`six>=1.5`; they do not need to be installed separately.

The exact environment used for verification was:

| Software | Version |
|---|---:|
| Python | 3.11.15 |
| NumPy | 2.4.6 |
| pandas | 3.0.3 |
| SciPy | 1.17.1 |
| python-dateutil | 2.9.0.post0 |
| six | 1.17.0 |
| pip | 24.0 |
| setuptools | 79.0.1 |

These exact pins are recorded in `requirements-tested.txt`. `pip check` was clean in the
verification environment.

### Hardware, memory, and storage

- **No non-standard hardware is required.** The demo, evaluator, and validity checks run
  on a CPU. A GPU is neither used nor detected.
- The real-data demo needs less than 200 MB of RAM and less than 1 MB of working storage.
- For the bundled benchmark tables, allow at least 2 GB of RAM for the largest
  pandas-based task and about 70 MB of disk after extraction. The compressed release
  archive is about 59 MB.

## 2. Installation guide

1. Extract the archive and enter its top-level directory:

   ```bash
   unzip EpiGym_v1.2_software.zip
   cd EpiGym_v1.2_software
   ```

2. Create and activate a virtual environment (POSIX shells):

   ```bash
   python3 -m venv .venv
   . .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install -r requirements.txt
   python -m pip check
   ```

   On Windows PowerShell, use `py -3 -m venv .venv`,
   `.venv\Scripts\Activate.ps1`, and then the same `pip install` commands.

3. No build step is required. The source files can be run directly with the environment's
   Python interpreter.

   ```bash
   python eval/evaluate.py --version
   # EpiGym evaluator 1.2
   ```

On a normal desktop computer with a broadband connection, creating the environment and
installing the dependencies typically takes **1-3 minutes**; extraction of the archive
normally takes less than one minute. A clean installation on the verification host took
17 seconds using cached binary wheels. Installation time is dominated by the local Python
package mirror and is not a property of EpiGym itself.

## 3. Demo

The demo uses a small, deterministic subset of the **real Tsuboyama stability benchmark**
and the corresponding real ProSST reference predictions. It exercises exactly the same T2
epistasis evaluator used for the full benchmark. The files contain 480 real double mutants
(120 from each of four domains), with measured score, additive prediction, measured
epistasis, and a ProSST prediction. The selection rule is recorded in
`demo/selection_manifest.json` and can be reproduced with
`python demo/make_demo_subset.py`; it does not inspect the resulting metric.

From the package root, run:

```bash
python demo/run_demo.py
```

Expected output is:

```text
{
  "task": "T2_epistasis_prediction",
  "protocol": "unsaturated subset [15,85]%",
  "median_abs_spearman": 0.3112,
  "mean_abs_spearman": 0.3064,
  "n_domains": 4,
  "n_expected": 275,
  "n_predicted": 275,
  "prediction_coverage": 1.0
}
DEMO PASS: output matches demo/expected_output.json
```

The demo normally completes in **2-5 seconds** on a normal desktop CPU and uses less than
200 MB of RAM. In the verification environment it took 1.9 seconds and used
135 MB peak resident memory. The numbers above are a software-path check on a small subset,
not a replacement for the full scientific benchmark or its leaderboard.

The equivalent command using the standalone evaluator is:

```bash
python eval/evaluate.py \
  --task T2_epistasis_prediction \
  --reference demo/demo_reference.csv \
  --pred demo/demo_predictions.csv
```

The optional `demo/make_demo_subset.py` script regenerates the two demo tables directly
from the bundled real benchmark and reference predictions:

```bash
python demo/make_demo_subset.py
```

## 4. Instructions for use

### Obtain the MAAD source data for T5

The MAAD antibody-neutralization source table is intentionally not included in this
archive because its upstream terms apply. Download it directly from the official source:

**https://www.raabmd.org/raab/download**

Prepare the downloaded records as an EpiGym-compatible reference table and supply that
table with `--reference`. The required T5 reference columns are `row_id`, `label`, and the
requested `split_*` columns; the frozen split names are listed in `tasks.json`. T1, T2, the
Demo, and the bundled validity checks do not require this download.

### Score a model

Run commands from the package root. A prediction file must contain one row per evaluated
variant and the key columns shown below.

```bash
# T1: high-order ProteinGym generalization
python eval/evaluate.py \
  --task T1_generalization \
  --split split_distance_controlled \
  --pred path/to/proteingym_predictions.csv

# T2: direct prediction of the measured epistatic residual
python eval/evaluate.py \
  --task T2_epistasis_prediction \
  --pred path/to/tsuboyama_predictions.csv

# T5: antibody-escape generalization
python eval/evaluate.py \
  --task T5_cross_assay_escape \
  --reference path/to/maad_escape.csv.gz \
  --pred path/to/maad_predictions.csv

# T5: evaluate selected splits only (optional)
python eval/evaluate.py \
  --task T5_cross_assay_escape \
  --reference path/to/maad_escape.csv.gz \
  --splits split_random,split_variant_holdout \
  --pred path/to/maad_predictions.csv
```

For T1 and T2, the evaluator reads the bundled reference table by default. T5 requires
the external reference table shown above. The evaluator prints a JSON result to standard
output; it does not train a model or download model weights.

### Prediction-file schemas

| Task | Required columns | Meaning |
|---|---|---|
| T1 | `assay_id`, `mutations`, `prediction` | predicted fitness; mutation strings use the bundled assay keys |
| T2 | `DMS_id`, `mutant`, `prediction` | predicted epistasis residual for a double mutant |
| T5 | `row_id`, `prediction` | escape score or probability |

Keys must match the corresponding bundled or user-supplied reference table. Prediction values must be
numeric and finite. By default, the evaluator requires a complete prediction file for the
requested evaluation subset and rejects missing rows, duplicate task keys, non-numeric
values, and infinite values. If a model cannot score every row, pass
`--allow-partial`; missing rows are then allowed but `n_expected`, `n_predicted`, and
`prediction_coverage` are included in the output and should be reported with the score.
Duplicate keys, non-numeric values, and infinite values are always rejected. Rows are
matched by the task key before scoring, and the frozen protocols are defined in
`tasks.json`.

```bash
python eval/evaluate.py \
  --task T2_epistasis_prediction \
  --allow-partial \
  --pred path/to/predictions_with_coverage_gaps.csv
```

### Use a small or external reference table

For a local test table with the same task columns, pass `--reference`:

```bash
python eval/evaluate.py \
  --task T2_epistasis_prediction \
  --reference path/to/reference.csv \
  --pred path/to/predictions.csv
```

For T1, the external table must also contain the requested split column (for example
`split_distance_controlled`). For T5, it must contain `label`, `row_id`, and any split
columns you want to score. The EpiGym MAAD reference schema defines 10 frozen `split_*`
columns. By default
T5 evaluates the six canonical splits reported in
`baselines/maad_escape_leaderboard.csv`; use `--splits` to select any other
comma-separated subset. The bundled T2 Demo illustrates the `--reference` option.

### Validation checks

To rerun the checks that the measured epistatic residual is structured signal rather than
measurement noise, run:

```bash
python eval/validity_checks.py
```

This command uses the bundled Tsuboyama and SKEMPI tables plus the frozen ProteinGym,
Tsuboyama, and MAAD reference summaries. It is separate from the small real-data Demo and
does not require the external MAAD source table.

## 5. Reproduction instructions (optional)

The evaluator self-test scores the three shipped T2 reference predictors through the same
code path used for the release leaderboard:

```bash
python eval/evaluate.py --selftest
```

Expected output covers all eight reference models:

```text
T2 ProSST    : median|rho|=0.3221  (n=50, coverage=100.0%)
T2 ESM2-650M : median|rho|=0.3195  (n=50, coverage=100.0%)
T2 PoET      : median|rho|=0.3157  (n=49, coverage=99.4%)
T2 ESM-C     : median|rho|=0.3077  (n=49, coverage=99.4%)
T2 EVmutation: median|rho|=0.2284  (n=47, coverage=88.9%)
T2 ESM1v     : median|rho|=0.2181  (n=50, coverage=100.0%)
T2 ESM2-3B   : median|rho|=0.2164  (n=50, coverage=100.0%)
T2 SPURS     : median|rho|=0.0949  (n=50, coverage=100.0%)
PASS: all eight models match baselines/tsuboyama_epistasis_leaderboard.csv
```

On the verification host, this self-test completed in 8.4 seconds on CPU and the full
validity check in 6.1 seconds. Both normally complete in under 15 seconds on a desktop.
Exact time varies with storage and the BLAS implementation. The package does not include
model weights or training code, so these commands reproduce evaluation and bundled
reference scores, not model inference.

## Contents

```text
EpiGym_v1.2_software/
├── README.md
├── requirements.txt
├── requirements-tested.txt
├── DATASHEET.md
├── VALIDITY.md
├── CHANGELOG.md
├── LICENSE.md
├── tasks.json
├── MANIFEST.json
├── CHECKSUMS.sha256
├── data/                 bundled benchmark tables (MAAD source table excluded)
├── baselines/            reference leaderboards and per-double predictions
├── eval/                 standalone evaluator and validity checks
└── demo/                 small real-data subset and runnable Demo
```

The complete task definitions and column-level provenance are in `tasks.json` and
`DATASHEET.md`. The frozen benchmark includes four assay types: ProteinGym DMS,
Tsuboyama stability, SKEMPI binding, and MAAD antibody escape. The MAAD source records are
user-supplied; the derived MAAD reference leaderboard remains bundled.

## Project repository

The canonical repository and issue tracker are at https://github.com/tubic/EpiGym. A
public release should use a versioned tag whose files and `CHECKSUMS.sha256` match this
release archive. Upstream data references are listed in `DATASHEET.md`.

## License and data notice

Derived EpiGym annotations, frozen splits, leaderboards, task cards, and evaluation code
are released under CC-BY-4.0.
