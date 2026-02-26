#!/usr/bin/env python3
"""
Count genes, exons, tRNA, and snoRNA features per chromosome from a GFF3 file.

Input:
  - GFF3 gzip file (data/saccharomyces_cerevisiae.gff.gz)
  - Chrom sizes file (data/chrom.sizes)

Output:
  - results/chr_feature_counts.tsv: Feature counts per chromosome
  - results/dropped_seqids.txt: seqids not found in chrom.sizes
"""

import gzip
from pathlib import Path
from collections import defaultdict
import pandas as pd

# Paths
repo_root = Path(__file__).parent.parent
gff_file = repo_root / "data" / "saccharomyces_cerevisiae.gff.gz"
chrom_sizes_file = repo_root / "data" / "chrom.sizes"
output_counts = repo_root / "results" / "chr_feature_counts.tsv"
output_dropped = repo_root / "results" / "dropped_seqids.txt"

# Read chrom.sizes
chrom_lengths = {}
with open(chrom_sizes_file) as f:
    for line in f:
        parts = line.strip().split()
        chrom, length = parts[0], int(parts[1])
        chrom_lengths[chrom] = length

valid_chroms = set(chrom_lengths.keys())

# Initialize counters
feature_counts = defaultdict(lambda: {
    "n_gene": 0,
    "n_exon_unique": set(),  # Use set to store unique (start, end, strand) tuples
    "n_tRNA": 0,
    "n_snoRNA": 0
})
dropped_seqids = set()
excluded_feature_lines = 0

# Read GFF3 and count features
with gzip.open(gff_file, "rt") as f:
    for line in f:
        # Skip comments and empty lines
        if line.startswith("#") or not line.strip():
            continue
        
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        
        seqid = parts[0]
        feature_type = parts[2]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        
        # Track dropped seqids and excluded lines
        if seqid not in valid_chroms:
            dropped_seqids.add(seqid)
            excluded_feature_lines += 1
            continue
        
        # Count by feature type
        if feature_type == "gene":
            feature_counts[seqid]["n_gene"] += 1
        elif feature_type == "exon":
            # Store unique (start, end, strand) tuples
            feature_counts[seqid]["n_exon_unique"].add((start, end, strand))
        elif feature_type == "tRNA":
            feature_counts[seqid]["n_tRNA"] += 1
        elif feature_type == "snoRNA":
            feature_counts[seqid]["n_snoRNA"] += 1

# Convert exon sets to counts and ensure all chromosomes are represented
for chrom in chrom_lengths:
    if chrom in feature_counts:
        feature_counts[chrom]["n_exon_unique"] = len(feature_counts[chrom]["n_exon_unique"])
    else:
        feature_counts[chrom] = {
            "n_gene": 0,
            "n_exon_unique": 0,
            "n_tRNA": 0,
            "n_snoRNA": 0
        }

# Save dropped seqids
with open(output_dropped, "w") as f:
    for seqid in sorted(dropped_seqids):
        f.write(f"{seqid}\n")

# Create output dataframe
rows = []
for chrom, length in sorted(chrom_lengths.items(), key=lambda x: x[1], reverse=True):
    counts = feature_counts[chrom]
    rows.append({
        "chrom": chrom,
        "chrom_length_bp": length,
        "n_gene": counts["n_gene"],
        "n_exon_unique": counts["n_exon_unique"],
        "n_tRNA": counts["n_tRNA"],
        "n_snoRNA": counts["n_snoRNA"]
    })

df = pd.DataFrame(rows)

# Calculate density columns (features per Mb)
df["gene_per_Mb"] = df["n_gene"] / (df["chrom_length_bp"] / 1e6)
df["exon_unique_per_Mb"] = df["n_exon_unique"] / (df["chrom_length_bp"] / 1e6)
df["tRNA_per_Mb"] = df["n_tRNA"] / (df["chrom_length_bp"] / 1e6)
df["snoRNA_per_Mb"] = df["n_snoRNA"] / (df["chrom_length_bp"] / 1e6)

# Round density columns to 4 decimal places
density_cols = ["gene_per_Mb", "exon_unique_per_Mb", "tRNA_per_Mb", "snoRNA_per_Mb"]
df[density_cols] = df[density_cols].round(4)

# Sort by gene_per_Mb descending
df = df.sort_values("gene_per_Mb", ascending=False).reset_index(drop=True)

# Save to TSV
df.to_csv(output_counts, sep="\t", index=False)

# Print statistics and results to console
print("=" * 80)
print("GFF3 Feature Counting Results")
print("=" * 80)
print(f"\nDropped seqids: {len(dropped_seqids)}")
print(f"Excluded feature lines: {excluded_feature_lines}")
print(f"\nTop 5 rows (sorted by gene_per_Mb descending):")
print(df.head(5).to_string(index=False))
print()
print(f"Full results saved to: {output_counts}")
print(f"Dropped seqids saved to: {output_dropped}")
