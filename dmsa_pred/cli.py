
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


# entry point
@group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Welcome to the dmsa-pred CLI!
    """
    pass


@cli.command(name="polyclonal-escape")
@option(
    "--alignment",
    required=True,
    type=Path(exists=True),
    help="",
)
@option(
    "--mut-escape-df",
    required=True,
    type=Path(exists=True),
    help="",
)
@option(
    "--activity-wt-df",
    required=True,
    type=Path(exists=True),
    help="",
)
@option(
    "--dms-wt-seq-id",
    required=True,
    type=str,
    help="",
)
@option(
    "--serum-label",
    required=True,
    type=str,
    help="",
)
@option(
    "--concentrations",
    required=True,
    type=str,
    help="",
)
@option(
    "--output",
    required=True,
    help="Path where the phip dataset will be dump'd to netCDF",
)
def polyclonal_escape_prediction(
    alignment,
    mut_escape_df,
    activity_wt_df,
    dms_wt_seq_id,
    serum_label,
    concentrations,
    output
):
    """
    """


    concentrations = [float(item) for item in concentrations_list.split(',')]
    mut_escape_df = pd.read_csv(mut_escape_df)

    # TODO remove these pre-emptively? or pass as parameter to script?
    sites_to_ignore = ["214a", "214b", "214c"]
    mut_escape_df = mut_escape_df[~mut_escape_df["site"].isin(sites_to_ignore)]
    mut_escape_df["escape"] = mut_escape_df["escape_median"]

    # Instantiate a Polyclonal object with betas and wildtype activity.
    model = polyclonal.Polyclonal(
        activity_wt_df=pd.read_csv(activity_wt_df),
        mut_escape_df=mut_escape_df,
        data_to_fit=None,
        alphabet=polyclonal.alphabets.AAS_WITHSTOP_WITHGAP,
    )

    # Mutation calling relative to the dms wildtype sequence.
    if alignment[-2:] == "xz":
        with lzma.open(alignment, "rt") as f:
            alignment = fasta_to_df(f)
    else:
        alignment = fasta_to_df(open(alignment, "r"))

    # TODO does this really even need to be in nextstrain tree?
    dms_wildtype = alignment.loc[dms_wt_seq, "seq"]
    
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
                f"prob_escape_{antibody}_c_{row.concentration}"
            ] = row.predicted_prob_escape

    write_json(ret_json, output)
