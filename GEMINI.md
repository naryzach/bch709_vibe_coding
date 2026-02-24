# Gemini Code Assist Configuration

Use this as a reference when interacting with Gemini in this repository.

## Conda environment
- The environment specification is in `environment.yml` at the project root.
- Environment name: `bch709_vibe_coding` using Python 3.11.
- Key Python packages: `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, `statsmodels`, `biopython`, `gffutils`, etc.
- There are also numerous R packages installed but they are not normally used by Python scripts.
- Create/activate with:
  ```bash
  conda env create -f environment.yml
  conda activate bch709_vibe_coding
  ```

## Repository layout
```
README.md
environment.yml

data/        # genomic and other input files (GFF, FASTA, chrom.sizes)
results/     # generated tables, plots, and other outputs
scripts/     # analytical code (Python, R, or notebooks)
```

Guidelines:
- Keep code within `scripts/`; do not sprinkle logic in notebooks or root.
- Treat `data/` as immutable; regenerate derived files into `results/`.

## Coding constraints
- Target Python 3.11; use modern syntax (f-strings, type hints).
- Adhere to PEP8 formatting; use linting if available.
- Avoid absolute paths; use `pathlib.Path` or `os.path` relative to repository.
- Guard execution with `if __name__ == '__main__'` when writing scripts.
- Use `pandas` dataframes for tabular data and `matplotlib`/`seaborn` for figures.

## Notes for Gemini
- The data directory contains sample yeast genomics data used in BCH709 homework and examples.
- The environment is comprehensive; rely on `environment.yml` for specifics if asked about available packages.
- When suggesting code, favor simple, educational examples and clearly comment their purpose.

This file is meant to be copy-pasted into Gemini's configuration or prompt to provide context about your workspace and constraints.