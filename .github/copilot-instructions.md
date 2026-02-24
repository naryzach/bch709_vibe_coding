# GitHub Copilot Instructions

This repository is used for BCH709 exercises and analyses.

## Conda environment
- Environment file: `environment.yml` at the root of the project.
- Name: `bch709_vibe_coding` (Python 3.11).
- Primary packages used in Python workflows: `pandas`, `numpy`, `matplotlib`, `seaborn`, `biopython`, `gffutils`, `scipy`, `statsmodels`, etc.
- Many additional R packages are also installed but Python code is the focus.

> **Tip**: always create/activate the conda env before running or writing code:
> ```sh
> conda env create -f environment.yml
> conda activate bch709_vibe_coding
> ```

## Project structure
```
README.md
environment.yml          # conda environment spec
data/                     # raw input files (compressed FASTA, GFF, chrom.sizes, etc.)
results/                  # output tables, figures, intermediate files
scripts/                  # user-written Python/R scripts or notebooks
```

- Scripts should live in `scripts/` and be designed to run from the repository root.
- Data files are read-only; never commit large generated results to source control.

## Code style
- Python 3.11 is the target; leverage f-strings, `pathlib.Path`, and type hints where practical.
- Follow PEP8/PEP257 formatting; use a formatter (black or `flake8`) if available.
- Docstrings are encouraged for public functions and scripts.
- Avoid hard‑coding absolute paths; build file paths relative to the repository root.
- Keep module imports at the top and avoid side effects; guard script execution with `if __name__ == "__main__":`.

## Architecture
- Very small and flat layout: data files live under `data/`, analysis code under `scripts/`, outputs go to `results/`.
- There is no build system or service layer; scripts should be runnable independently from the repo root.
- Data flows are simple: read from `data/`, process in Python/R, write tabular results or figures to `results/`.

## Build and test
- No compilation or build step. To run code use `python scripts/<name>.py` after activating the conda environment.
- The environment is built with `conda env create -f environment.yml`.
- There is currently no automated test suite; if tests are added they should live alongside scripts and use `pytest`.

## Project conventions
- All new Python or R code belongs in `scripts/`; avoid mixing logic into the root.
- Treat data in `data/` as read-only. Derived files, intermediate tables and plots should be written to `results/` and are typically ignored by git.
- Solutions should be simple and explicit—this is an educational repo, not production software.
- Update `environment.yml` when adding new package dependencies.

## Integration points
- External data is fetched manually (see download commands in existing terminal history) and stored in `data/`.
- There are no APIs, databases, or external services; analyses operate on local files only.

## Security
- No credentials, tokens, or secrets are present anywhere in the repository.
- The data is public yeast genomics data; there are no privacy concerns.

This document helps Copilot know how to suggest context‑appropriate completions. Refer to the environment.yml for a full list of packages and versions.