import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def main():
    parser = argparse.ArgumentParser(description="Plot KDE density of overlap vs log10(HMM e-value)")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/domain_comparison_results_human.tsv"),
        help="Path to TSV input file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("density_human.png"),
        help="Path to output PNG",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    df['log_hmm_evalue'] = np.log10(df['hmm_evalue'].astype(float) + 1e-300)

    plt.figure(figsize=(10, 6))
    sns.kdeplot(
        x=df['log_hmm_evalue'],
        y=df['overlap_percentage'],
        fill=True,
        cmap='viridis',
        bw_adjust=0.5,
        levels=100,
    )

    plt.xlabel('Log10(Evalue)')
    plt.ylabel('Overlap Percentage')
    plt.title('Density Plot of Overlap vs Log10(Evalue)')

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
