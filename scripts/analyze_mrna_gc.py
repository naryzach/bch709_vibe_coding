"""
Script to analyze GC content distribution of yeast mRNA sequences.

This script:
1. Downloads mRNA FASTA data from UCSC if missing
2. Parses multi-line FASTA records directly from gzip compression
3. Computes sequence length and GC content for each mRNA
4. Generates summary statistics and visualizations
5. Outputs metrics table and distribution plot

Input specification:
- File: data/mrna.fa.gz (FASTA, gzip-compressed)
  - Download from:
    https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/mrna.fa.gz
  - Header format:
    >ACCESSION additional metadata
  - Accession ID = first token after ">"
  - Sequences may span multiple lines

Analysis conditions:
- Read mrna.fa.gz directly using gzip (no manual decompression)
- Parse multi-line FASTA records correctly

For each sequence:
- Extract accession ID
- Compute sequence length = count(A, T, G, C, N)
- Compute GC content:
  GC = (count(G) + count(C)) / (count(A) + count(T) + count(G) + count(C))
- Ignore ambiguous bases (N) in GC denominator
- Remove sequences with:
  - Zero valid (A/T/G/C) bases
  - All NA values

Output 1 – Table:
- Filename: results/mrna_metrics.tsv
- Columns:
  - accession
  - length
  - gc_content
- Decimal places: gc_content = 4
- Sorting: gc_content descending
- Filter: none (include all valid sequences)

Output 2 – Plot:
- Filename: results/gc_content_distribution.png
- Size: 1600 × 900 px, dpi: 200
- Colors:
  - Histogram: default matplotlib/seaborn color
  - Density curve: contrasting solid color
  - Mean/median lines: black dashed
- Axes / labels / legend:
  - X-axis: GC content (0–1)
  - Y-axis: Density
  - Title: "GC Content Distribution of Yeast mRNA"
  - Label axes clearly
  - Include legend for mean/median
- Plot elements:
  - Histogram of GC values
  - Kernel density overlay
  - Vertical dashed lines for mean and median
  - Caption showing:
    n, mean, median, sd (4 decimals)

Output 3 – QC:
- Dropped items file: none
- Print to console:
  - Total sequences read
  - Valid sequences analyzed
  - Mean, median, SD of GC content

Additional requirements:
- Automatically download mrna.fa.gz if missing
- Create directories if missing:
  data/
  results/
- Use main() function
- Handle I/O errors gracefully
- No user input required
- Script must run with:
  python analyze_mrna_gc.py
- Ensure reproducible output
"""

import gzip
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import statistics

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

try:
    import urllib.request
except ImportError:
    import urllib2 as urllib


def ensure_directories() -> None:
    """Create data/ and results/ directories if missing."""
    Path("data").mkdir(exist_ok=True)
    Path("results").mkdir(exist_ok=True)


def download_mrna_fasta(url: str, filepath: str) -> None:
    """Download mRNA FASTA file if not present."""
    if Path(filepath).exists():
        print(f"✓ {filepath} already exists")
        return
    
    print(f"⬇ Downloading mRNA FASTA from {url}...")
    try:
        urllib.request.urlretrieve(url, filepath)
        print(f"✓ Downloaded to {filepath}")
    except Exception as e:
        print(f"✗ Error downloading file: {e}", file=sys.stderr)
        raise


def parse_fasta_gz(filepath: str) -> Dict[str, str]:
    """
    Parse multi-line FASTA records from gzip-compressed file.
    
    Args:
        filepath: Path to gzip-compressed FASTA file
        
    Returns:
        Dictionary mapping accession ID to sequence
    """
    sequences = {}
    current_accession = None
    current_seq = []
    
    with gzip.open(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\n')
            
            if not line:
                continue
            
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_accession is not None:
                    sequences[current_accession] = ''.join(current_seq).upper()
                
                # Extract accession ID (first token after '>')
                header = line[1:].split()
                current_accession = header[0] if header else f"seq_{len(sequences)}"
                current_seq = []
            else:
                # Accumulate sequence lines
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_accession is not None:
            sequences[current_accession] = ''.join(current_seq).upper()
    
    return sequences


def compute_gc_content(seq: str) -> Tuple[float, int]:
    """
    Compute GC content and sequence length.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        Tuple of (gc_content, length)
        - gc_content: fraction between 0 and 1 (or NaN if no valid bases)
        - length: total count of A, T, G, C, N
    """
    # Count bases (case-insensitive, uppercase already)
    a_count = seq.count('A')
    t_count = seq.count('T')
    g_count = seq.count('G')
    c_count = seq.count('C')
    n_count = seq.count('N')
    
    # Total length (A + T + G + C + N)
    total_length = a_count + t_count + g_count + c_count + n_count
    
    # Valid bases for GC calculation (A + T + G + C, exclude N)
    valid_bases = a_count + t_count + g_count + c_count
    
    # GC content
    if valid_bases > 0:
        gc_content = (g_count + c_count) / valid_bases
    else:
        gc_content = np.nan
    
    return gc_content, total_length


def analyze_mrna_metrics(filepath: str) -> pd.DataFrame:
    """
    Analyze mRNA sequences and compute GC content metrics.
    
    Args:
        filepath: Path to gzip-compressed FASTA file
        
    Returns:
        DataFrame with columns: accession, length, gc_content
    """
    print(f"\n📖 Parsing {filepath}...")
    sequences = parse_fasta_gz(filepath)
    print(f"✓ Read {len(sequences)} sequences")
    
    records = []
    for accession, seq in sequences.items():
        gc_content, length = compute_gc_content(seq)
        
        # Skip sequences with zero valid bases or all NaN
        if length == 0 or np.isnan(gc_content):
            continue
        
        records.append({
            'accession': accession,
            'length': length,
            'gc_content': gc_content
        })
    
    df = pd.DataFrame(records)
    print(f"✓ Valid sequences: {len(df)}")
    print(f"✓ Dropped: {len(sequences) - len(df)} sequences")
    
    return df


def compute_statistics(gc_values: np.ndarray) -> Dict[str, float]:
    """Compute summary statistics for GC content."""
    return {
        'n': len(gc_values),
        'mean': np.mean(gc_values),
        'median': np.median(gc_values),
        'sd': np.std(gc_values, ddof=1)  # Sample SD
    }


def create_distribution_plot(gc_values: np.ndarray, stats: Dict[str, float], 
                            output_path: str) -> None:
    """
    Create histogram with density overlay and statistics.
    
    Args:
        gc_values: Array of GC content values
        stats: Dictionary with n, mean, median, sd
        output_path: Path to save PNG file
    """
    # Set figure size and DPI to achieve 1600x900 px at 200 DPI
    # 1600 px / 200 dpi = 8 inches
    # 900 px / 200 dpi = 4.5 inches
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=200)
    
    # Histogram with density normalization
    ax.hist(gc_values, bins=50, density=True, alpha=0.7, 
            edgecolor='black', linewidth=0.5, label='Histogram')
    
    # Kernel Density Estimate overlay
    from scipy import stats as sp_stats
    kde = sp_stats.gaussian_kde(gc_values)
    x_range = np.linspace(0, 1, 500)
    ax.plot(x_range, kde(x_range), color='red', linewidth=2, label='Density')
    
    # Mean and median lines
    mean_val = stats['mean']
    median_val = stats['median']
    ax.axvline(mean_val, color='black', linestyle='--', linewidth=1.5, 
               label=f"Mean ({mean_val:.4f})")
    ax.axvline(median_val, color='black', linestyle='--', linewidth=1.5, 
               label=f"Median ({median_val:.4f})")
    
    # Labels and title
    ax.set_xlabel('GC Content', fontsize=12, fontweight='bold')
    ax.set_ylabel('Density', fontsize=12, fontweight='bold')
    ax.set_title('GC Content Distribution of Yeast mRNA', fontsize=14, fontweight='bold')
    ax.set_xlim(0, 1)
    
    # Caption with statistics
    caption = (f"n = {stats['n']}, "
               f"mean = {stats['mean']:.4f}, "
               f"median = {stats['median']:.4f}, "
               f"SD = {stats['sd']:.4f}")
    ax.text(0.98, 0.97, caption, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Legend
    ax.legend(loc='upper left', fontsize=10)
    
    # Grid for readability
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"✓ Saved plot to {output_path}")
    plt.close()


def main() -> None:
    """Main analysis pipeline."""
    print("=" * 70)
    print("Yeast mRNA GC Content Analysis")
    print("=" * 70)
    
    try:
        # Setup
        ensure_directories()
        
        # Download if needed
        url = "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/mrna.fa.gz"
        fasta_file = "data/mrna.fa.gz"
        download_mrna_fasta(url, fasta_file)
        
        # Analyze sequences
        df = analyze_mrna_metrics(fasta_file)
        
        # Sort by GC content descending
        df = df.sort_values('gc_content', ascending=False).reset_index(drop=True)
        
        # Round GC content to 4 decimals
        df['gc_content'] = df['gc_content'].round(4)
        
        # Save metrics table
        output_tsv = "results/mrna_metrics.tsv"
        df.to_csv(output_tsv, sep='\t', index=False)
        print(f"✓ Saved metrics to {output_tsv}")
        
        # Compute statistics
        gc_values = df['gc_content'].values
        stats = compute_statistics(gc_values)
        
        # Print QC summary
        print("\n" + "=" * 70)
        print("Quality Control Summary")
        print("=" * 70)
        print(f"Total sequences read:       {len(df)}")
        print(f"Valid sequences analyzed:   {len(df)}")
        print(f"Mean GC content:            {stats['mean']:.4f}")
        print(f"Median GC content:          {stats['median']:.4f}")
        print(f"Std Dev GC content:         {stats['sd']:.4f}")
        print("=" * 70)
        
        # Create visualization
        output_plot = "results/gc_content_distribution.png"
        create_distribution_plot(gc_values, stats, output_plot)
        
        # Biological interpretation
        print("\n" + "=" * 70)
        print("Biological Interpretation")
        print("=" * 70)
        interpretation = """
TWO biological reasons for the observed GC-content distribution in yeast mRNA:

1. **Codon Usage Bias and Translational Efficiency**
   Yeast genes exhibit strong codon preference driven by tRNA availability and 
   ribosomal efficiency. Genes encoding abundant proteins (ribosomal proteins, 
   metabolic enzymes) preferentially use codons matching abundant tRNAs, which 
   correlates with specific nucleotide compositions. GC-rich codons may be 
   selected or avoided based on optimal translation speed and accuracy for 
   highly-expressed genes, creating distinct GC content clusters. The bimodal 
   or skewed distribution observed may reflect different gene expression level 
   classes (highly-expressed vs. lowly-expressed), each with optimized codon 
   usage patterns.

2. **mRNA Stability and Selection Pressure on Thermodynamic Properties**
   GC-rich mRNA sequences form more stable secondary structures (G-C base pairs 
   have 3 hydrogen bonds vs. 2 for A-T), affecting mRNA half-life and protein 
   synthesis rates. Genes requiring stable transcripts (e.g., those encoding 
   structural proteins or long-lasting regulatory proteins) may be selected for 
   higher GC content to resist degradation in cellular conditions. Conversely, 
   genes requiring rapid turnover (e.g., cell cycle regulators, stress response 
   genes) may evolve lower GC content for shorter half-lives. This selective 
   pressure, combined with the constraint that GC content must remain compatible 
   with functional codon usage patterns, creates the observed distribution 
   reflecting a balance between translation efficiency and mRNA stability needs.
        """
        print(interpretation)
        
        print("\n✓ Analysis complete!")
        
    except Exception as e:
        print(f"\n✗ Error during analysis: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
