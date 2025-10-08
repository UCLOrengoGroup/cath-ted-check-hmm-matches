import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLOR_MAPPING = {'TED_better_hmm_evalue': 'red', 'hmm_better_evalue': 'black'}


def elbow_point(x: np.ndarray, y: np.ndarray) -> int:
    coords = np.vstack((x, y)).T
    origin = coords[0]
    vec = coords[-1] - origin
    norm = np.sqrt(np.sum(vec ** 2))
    if norm == 0:
        return 0
    vec = vec / norm
    vec_from_origin = coords - origin
    scalar_product = np.sum(vec_from_origin * vec, axis=1)
    proj = np.outer(scalar_product, vec)
    dist = np.sqrt(np.sum((vec_from_origin - proj) ** 2, axis=1))
    return int(np.argmax(dist))


def main():
    parser = argparse.ArgumentParser(description="Scatterplot and derived plots for TED e-values vs overlap")
    parser.add_argument("--input", type=Path, default=Path("ted-evalue.txt"), help="Whitespace-delimited input file")
    parser.add_argument("--outdir", type=Path, default=Path("."), help="Directory to write outputs")
    parser.add_argument("--ymax", type=float, default=80.0, help="Max y-axis for filtered plots")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(args.input, delim_whitespace=True)

    filtered = data[data['overlap_percentage'] < args.ymax]

    # Scatter filtered
    plt.figure(figsize=(12, 8))
    plt.scatter(
        filtered['e-value'],
        filtered['overlap_percentage'],
        c=filtered['source'].map(COLOR_MAPPING),
        s=80,
        alpha=0.7,
        edgecolors='black',
    )
    plt.xscale('log')
    plt.xlabel('E-value', fontsize=14)
    plt.ylabel('Overlap (%)', fontsize=14)
    plt.title(f'E-value vs Overlap Scatterplot (Overlap < {args.ymax}%)', fontsize=16)
    plt.ylim(max(0, filtered['overlap_percentage'].min() - 5), args.ymax)
    plt.tight_layout()
    plt.savefig(args.outdir / 'evalue_vs_overlap_scatterplot_filtered.png', dpi=300)
    plt.close()

    # Histogram
    plt.figure(figsize=(12, 6))
    plt.hist(filtered['e-value'], bins=50, edgecolor='black')
    plt.xscale('log')
    plt.xlabel('E-value', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title('Histogram of E-values (filtered by overlap)', fontsize=16)
    median_evalue = np.median(filtered['e-value'])
    mean_evalue = np.mean(filtered['e-value'])
    plt.axvline(median_evalue, color='red', linestyle='dashed', linewidth=2, label=f'Median: {median_evalue:.2e}')
    plt.axvline(mean_evalue, color='green', linestyle='dashed', linewidth=2, label=f'Mean: {mean_evalue:.2e}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.outdir / 'evalue_histogram_filtered.png', dpi=300)
    plt.close()

    # CDF
    plt.figure(figsize=(12, 6))
    sorted_evalues = np.sort(filtered['e-value'])
    yvals = np.arange(len(sorted_evalues)) / float(max(len(sorted_evalues) - 1, 1))
    plt.plot(sorted_evalues, yvals)
    plt.xscale('log')
    plt.xlabel('E-value', fontsize=14)
    plt.ylabel('Cumulative Probability', fontsize=14)
    plt.title('CDF of E-values (filtered by overlap)', fontsize=16)
    plt.tight_layout()
    plt.savefig(args.outdir / 'evalue_cdf_filtered.png', dpi=300)
    plt.close()

    # Elbow
    x = np.log10(sorted_evalues)
    y = np.arange(1, len(sorted_evalues) + 1) / len(sorted_evalues)
    idx = elbow_point(x, y) if len(sorted_evalues) > 1 else 0
    elbow_evalue = sorted_evalues[idx] if len(sorted_evalues) else float('nan')

    # Scatter with cutoffs
    plt.figure(figsize=(12, 8))
    plt.scatter(
        filtered['e-value'],
        filtered['overlap_percentage'],
        c=filtered['source'].map(COLOR_MAPPING),
        s=80,
        alpha=0.7,
        edgecolors='black',
    )
    plt.xscale('log')
    plt.xlabel('E-value', fontsize=14)
    plt.ylabel('Overlap (%)', fontsize=14)
    plt.title('E-value vs Overlap with Potential Cutoffs', fontsize=16)
    plt.ylim(max(0, filtered['overlap_percentage'].min() - 5), args.ymax)
    if not np.isnan(elbow_evalue):
        plt.axvline(np.median(filtered['e-value']), color='red', linestyle='dashed', linewidth=2, label=f'Median: {median_evalue:.2e}')
        plt.axvline(np.mean(filtered['e-value']), color='green', linestyle='dashed', linewidth=2, label=f'Mean: {mean_evalue:.2e}')
        plt.axvline(elbow_evalue, color='blue', linestyle='dashed', linewidth=2, label=f'Elbow: {elbow_evalue:.2e}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.outdir / 'ted-evalue_vs_overlap_scatterplot_with_cutoffs_filtered.png', dpi=300)
    plt.close()

    print("Potential E-value cutoffs (filtered):")
    print(f"Median: {median_evalue:.2e}")
    print(f"Mean: {mean_evalue:.2e}")
    if not np.isnan(elbow_evalue):
        print(f"Elbow point: {elbow_evalue:.2e}")


if __name__ == "__main__":
    main()
