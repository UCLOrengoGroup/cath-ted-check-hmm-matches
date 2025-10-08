#!/usr/bin/env python3

import pandas as pd

# We're going to add annotation columns to this results file ...
domain_path = "results/domain_comparison_results_human.processed.tsv"

# ... and create this file
results_path = "results/domain_comparison_results_human.processed.annotated.tsv"

# using these files
anno_colnames = "ted_id	md5_domain	consensus_level	chopping	nres_domain	num_segments	plddt	num_helix_strand_turn	num_helix	num_strand	num_helix_strand	num_turn	proteome_id	cath_label	cath_level	cath_method	packing_density	norm_rg	tax_common	tax_scientific	tax_lineage".split(
    "\t"
)
anno_path = "ted_100_324m.domain_summary.cath.globularity.taxid.human.tsv"


def run():

    # read in original results file
    domain_df = pd.read_csv(domain_path, sep="\t")
    domain_df["ted_id"] = domain_df["foldseek_domain_id"]
    domain_df.set_index("ted_id", inplace=True)

    # get plddt / packing
    anno_df = pd.read_csv(anno_path, names=anno_colnames, sep="\t")
    anno_df.set_index("ted_id", inplace=True)

    # add plddt / packing
    output_df = domain_df.join(anno_df[["plddt", "packing_density"]], how="inner")

    output_df.to_csv(results_path, sep="\t", index=False)

    print(output_df.head())


if __name__ == "__main__":
    run()
