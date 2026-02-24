# ChatGPT / Codex Context Document

This markdown file summarizes the project environment, structure and coding conventions for use when working with ChatGPT or Codex on the repository.

## Environment overview
- Conda environment file: `environment.yml` (root of repo).
- Environment name: `bch709_vibe_coding` (Python 3.11).
- Contains numerous scientific packages: e.g. `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, `statsmodels`, `biopython`, `gffutils` plus many others.
- R packages are also installed but the focus is Python-based analyses.
- To use the environment:
  ```bash
  conda env create -f environment.yml
  conda activate bch709_vibe_coding
  ```

## Project layout
```
README.md

environment.yml          # conda environment specification

data/                    # raw inputs (FASTAs, GFFs, chrom sizes, expression tables)
results/                 # output tables, plots, and derived files
scripts/                 # Python/R scripts or Jupyter notebooks containing analysis code
```
- All new code should go under `scripts/` and be runnable from the repository root.
- Data in `data/` should remain unchanged; generated results go to `results/`.

## Coding guidelines
- Target Python 3.11; use modern features like f-strings and type hints.
- Follow PEP8/PEP257 style, add docstrings and comments where helpful.
- Use relative paths for file access (`pathlib.Path` or `os.path.join`).
- Structure scripts with a `main()` and guard `if __name__ == '__main__':`.
- Prefer `pandas` for tabular manipulation and `matplotlib`/`seaborn` for plotting.

## Additional notes for Codex
- The repository is intended for teaching; solutions should be clear and explicit.
- Avoid introducing dependencies not listed in `environment.yml` unless they are trivial and can be added by updating the environment.
- Refer back to `environment.yml` to answer questions about installed package versions.

This document is for copy‑paste into prompts or notes when collaborating with Codex or similar AI assistants.