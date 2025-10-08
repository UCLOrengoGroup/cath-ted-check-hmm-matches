#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

base_dir = Path().parent.parent.resolve()

data_file = base_dir / "data/domain_comparison_results_human.tsv"
old_processed_file = base_dir / "results/domain_comparison_results_human.processed.tsv"
new_processed_file = (
    base_dir / "results/domain_comparison_results_human.processed.tsv.new"
)


def run():
    # Load the data
    df_data = pd.read_csv(str(data_file), sep="\t")
    df_data["ted_id"] = df_data["foldseek_domain_id"]
    df_data.set_index("ted_id", inplace=True)

    # Load the old processed data
    df_old = pd.read_csv(str(old_processed_file), sep="\t")
    df_old["ted_id"] = df_old["foldseek_domain_id"]
    df_old.set_index("ted_id", inplace=True)

    # Join the new column to the existing data
    df_new = df_old.join(df_data[["hmm_id"]], how="inner")

    # Save the new processed data
    df_new.to_csv(str(new_processed_file), sep="\t", index=False)

    # Do the Indy swap on success
    old_processed_file.rename(old_processed_file.with_suffix(".bak"))
    new_processed_file.rename(old_processed_file)


if __name__ == "__main__":
    run()
