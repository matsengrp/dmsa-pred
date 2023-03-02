
"""
@File: dmsa-pred_cli.py

@Author: Jared Galloway

Command line interface (CLI) for dmsa-pred.
"""

# built-in
import gzip
import glob
import sys
import lzma
from collections import defaultdict
import json
import natsort
# import argparse

# dependencies
import pandas as pd
from click import Path, group, option, argument
import click

# existing deps
from augur.utils import write_json
from Bio import SeqIO
from Bio.Seq import Seq

# new deps
import polyclonal

# shared utilities
# TODO Bio.Align may be faster?
def fasta_to_df(fasta_file):
    """simply convert a fasta to dataframe"""

    ids, seqs = [], []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):  # (generator)
        ids.append(seq_record.id)
        seqs.append(str(seq_record.seq))
    return pd.DataFrame({"strain": ids, "seq": seqs}).set_index("strain")


def mutations(naive_aa, aa, allowed_subs):
    """Amino acid substitutions between two sequences, in IMGT coordinates."""

    assert len(naive_aa) == len(aa)
    return " ".join(
        [
            f"{aa1}{pos+1}{aa2}"
            for pos, (aa1, aa2) in enumerate(zip(naive_aa, aa))
            if aa1 != aa2 and f"{aa1}{pos+1}{aa2}" in allowed_subs
        ]
    )

# define all common options
def group_options(*options):
    """
    custom decorator for shared options
    """
    def wrapper(function):
        for option in reversed(options):
            function = option(function)
        return function
    return wrapper

alignment = option(
    "--alignment",
    required=True,
    type=Path(exists=True),
    help="",
)
mut_effects_df = option(
    "--mut-effects-df",
    required=True,
    type=Path(exists=True),
    help="",
)
dms_wt_seq_id = option(
    "--dms-wt-seq-id",
    required=True,
    type=str,
    help="",
)
experiment_label = option(
    "--experiment-label",
    required=True,
    type=str,
    help="",
)
output = option(
    "--output",
    required=True,
    help="Path where the phip dataset will be dump'd to netCDF",
)

# entry point
@group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Welcome to the dmsa-pred CLI!
    """
    pass

@cli.command(name="polyclonal-escape")
@option(
    "--activity-wt-df",
    required=True,
    type=Path(exists=True),
    help="",
)
@option(
    "--concentrations",
    required=True,
    type=str,
    help="",
)
@option(
    "--escape-column",
    required=False,
    default="escape",
    type=str,
    help="",
)
@group_options(alignment, mut_effects_df, dms_wt_seq_id, experiment_label, output)
def polyclonal_escape_prediction(
    activity_wt_df,
    concentrations,
    escape_column,
    alignment,
    mut_effects_df,
    dms_wt_seq_id,
    experiment_label,
    output
):
    """
    """

    concentrations = [float(item) for item in concentrations.split(',')]
    mut_effects_df = pd.read_csv(mut_effects_df).rename({escape_column:"escape"}, axis=1)

    # TODO remove these pre-emptively? or pass as parameter to script?
    #sites_to_ignore = ["214a", "214b", "214c"]
    #mut_effects_df = mut_effects_df[~mut_effects_df["site"].isin(sites_to_ignore)]
    #mut_effects_df["escape"] = mut_effects_df["escape_median"]

    # Instantiate a Polyclonal object with betas and wildtype activity.
    model = polyclonal.Polyclonal(
        activity_wt_df=pd.read_csv(activity_wt_df),
        mut_escape_df=mut_effects_df,
        data_to_fit=None,
        sites=tuple(natsort.natsorted(mut_effects_df.site, alg=natsort.ns.SIGNED)),
        alphabet=polyclonal.alphabets.AAS_WITHSTOP_WITHGAP,
    )

    # Mutation calling relative to the dms wildtype sequence.
    if alignment[-2:] == "xz":
        with lzma.open(alignment, "rt") as f:
            alignment = fasta_to_df(f)
    else:
        alignment = fasta_to_df(open(alignment, "r"))

    # TODO does this really even need to be in nextstrain tree?
    dms_wildtype = alignment.loc[dms_wt_seq_id, "seq"]
    
    # TODO N jobs? pandarallel apply()
    alignment["aa_substitutions"] = alignment.seq.apply(
        lambda aligned_seq: mutations(dms_wildtype, aligned_seq, set(model.mutations))
    )
    alignment.reset_index(inplace=True)

    escape_probs = (
        model.prob_escape(
            variants_df=alignment, concentrations=concentrations
        )
        .drop("seq", axis=1)
        .reset_index()
    )

    ret_json = {"generated_by": {"program": "polyclonal"}, "nodes": defaultdict(dict)}
    for strain, strain_df in escape_probs.groupby("strain"):
        for idx, row in strain_df.iterrows():
            ret_json["nodes"][strain][
                f"prob_escape_{experiment_label}_c_{row.concentration}"
            ] = row.predicted_prob_escape

    write_json(ret_json, output)

@cli.command(name="escape-fraction")
@group_options(alignment, mut_effects_df, dms_wt_seq_id, experiment_label, output)
def escape_fraction_prediction(
    alignment,
    mut_effects_df,
    dms_wt_seq_id,
    experiment_label,
    output
):
    """
    """

    mut_effects_df = pd.read_csv(mut_effects_df)
    mut_effects_df = mut_effects_df.assign(
        non_escape_frac=(1-mut_effects_df["mut_escape_frac_epistasis_model"])
    )

    if alignment[-2:] == "xz":
        with lzma.open(alignment, "rt") as f:
            alignment = fasta_to_df(f)
    else:
        alignment = fasta_to_df(open(alignment, "r"))

    dms_wildtype = alignment.loc[dms_wt_seq_id, "seq"]
    
    alignment["aa_substitutions"] = alignment.seq.apply(
        lambda aligned_seq: mutations(dms_wildtype, aligned_seq, set(mut_effects_df.aa_substitution))
    )
    alignment.reset_index(inplace=True)

    def compute_variant_escape_score(
        aa_subs,
        mut_effect_col = "non_escape_frac"
    ):
        data = mut_effects_df[mut_effects_df['aa_substitution'].isin(aa_subs.split())]
        return 1-data[mut_effect_col].prod()

    alignment["variant_escape_score"] = alignment.aa_substitutions.apply(
        lambda aa_subs: compute_variant_escape_score(aa_subs)
    )

    ret_json = {"generated_by": {"program": "custom"}, "nodes": defaultdict(dict)}
    for idx, row in alignment.iterrows():
        ret_json["nodes"][row.strain][
            f"{experiment_label}_escape_fraction"
        ] = row.variant_escape_score

    write_json(ret_json, output)

@cli.command(name="additive-phenotype")
@option(
    "--phenotype-column",
    required=False,
    default="escape",
    type=str,
    help="",
)
@group_options(alignment, mut_effects_df, dms_wt_seq_id, experiment_label, output)
def additive_phenotype_prediction(
    phenotype_column,
    alignment,
    mut_effects_df,
    dms_wt_seq_id,
    experiment_label,
    output
):
    """
    """

    mut_effects_df = pd.read_csv(mut_effects_df)

    if alignment[-2:] == "xz":
        with lzma.open(alignment, "rt") as f:
            alignment = fasta_to_df(f)
    else:
        alignment = fasta_to_df(open(alignment, "r"))

    dms_wildtype = alignment.loc[dms_wt_seq_id, "seq"]
    
    alignment["aa_substitutions"] = alignment.seq.apply(
        lambda aligned_seq: mutations(dms_wildtype, aligned_seq, set(mut_effects_df.aa_substitution))
    )
    alignment.reset_index(inplace=True)

    def compute_variant_additive_phenotype(
        aa_subs
    ):
        data = mut_effects_df[mut_effects_df['aa_substitution'].isin(aa_subs.split())]
        return data[phenotype_column].sum()

    alignment["variant_escape_score"] = alignment.aa_substitutions.apply(
        lambda aa_subs: compute_variant_additive_phenotype(aa_subs)
    )

    ret_json = {"generated_by": {"program": "custom"}, "nodes": defaultdict(dict)}
    for idx, row in alignment.iterrows():
        ret_json["nodes"][row.strain][
            f"{experiment_label}_additive_phenotype"
        ] = row.variant_escape_score

    write_json(ret_json, output)

if __name__ == '__main__':
    cli()
