"""
@File: dmsa-pred_cli.py

@Author: Jared Galloway and Hugh Haddox

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
def parse_fasta_entries(open_fasta_file):
    """Parse entries in a FASTA file, storing each in rows of a dataframe"""
    
    ids, seqs = [], []
    for seq_record in SeqIO.parse(open_fasta_file, "fasta"):
        ids.append(seq_record.id)
        seqs.append(str(seq_record.seq))
    return pd.DataFrame({"strain": ids, "seq": seqs}).set_index("strain")

def fasta_to_df(fasta_file):
    """Open a FASTA file, then read in each entry and store in a dataframe"""
    
    if fasta_file[-2:] == "xz":
        with lzma.open(fasta_file, "rt") as f:
            return parse_fasta_entries(f)
    else:
        with open(fasta_file) as f:
            return parse_fasta_entries(f)

def get_mutations(seq1, seq2, allowed_mutations):
    """Make a list of amino-acid mutations between two sequences"""

    assert len(seq1) == len(seq2)
    return " ".join(
        [
            f"{aa1}{pos+1}{aa2}"
            for pos, (aa1, aa2) in enumerate(zip(seq1, seq2))
            if aa1 != aa2 and f"{aa1}{pos+1}{aa2}" in allowed_mutations
        ]
    )

def write_output_json(alignment_df, phenotype_col, output_json):
    """
    Write an output JSON for use in a Nextstrain workflow, storing DMS
    predictions for each sequence in the alignment.
    
    Loops over entries in the dataframe and records predicted phenotypes,
    grouping by strain to account for dataframes with multiple
    predictions per strain (e.g., different polyclonal concentrations)
    """
    
    ret_json = {
        "generated_by": {"program": "dmsa-pred"},
        "nodes": defaultdict(dict)
    }
    for strain, strain_df in alignment_df.groupby("strain"):
        #ret_json["nodes"][strain]["DMS"] = 
        for idx, row in strain_df.iterrows():
            ret_json["nodes"][strain][row['json_label']] = row[phenotype_col]
    write_json(ret_json, output_json)

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
    help="a FASTA file with a multiple-sequence alignment of input sequences",
)
dms_wt_seq_id = option(
    "--dms-wt-seq-id",
    required=True,
    type=str,
    help="the name of the DMS wildtype sequence from the input FASTA file",
)
mut_effects_df = option(
    "--mut-effects-df",
    required=True,
    type=Path(exists=True),
    help="a dataframe giving the effects of mutations",
)
mut_effect_col = option(
    "--mut-effect-col",
    required=True,
    type=str,
    help="a column from mut_effects_df that gives the effect of the mutation",
)
# TODO remove this
site_col = option(
    "--site-col",
    required=True,
    type=str,
    help="a column from mut_effects_df that gives the site of the mutation",
)
experiment_label = option(
    "--experiment-label",
    required=True,
    type=str,
    help="",
)
output_json = option(
    "--output-json",
    required=True,
    help="an output JSON file for use in a Nextstrain workflow",
)
output_df = option(
    "--output-df",
    required=True,
    help="an output dataframe with mutations and DMS predictions for each sequence in the input alignment",
)

# entry point
@group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Welcome to the dmsa-pred CLI!
    """
    pass


#@cli.command(name="call-mutations")
#@group_options(alignment, dms_wt_seq_id, output_json, output_df)
#def call_mutations(
#    alignment,
#    dms_wt_seq_id,
#    output_json,
#    output_df
#):
#    """
#    Read in the input alignment FASTA. Then, for each sequence in the
#    # alignment, make a list of all mutations relative to the DMS WT sequence
#    """
#
#    alignment_df = fasta_to_df(alignment)
#    dms_wt_seq = alignment_df.loc[dms_wt_seq_id, "seq"]
#    alignment_df.reset_index(inplace=True)
#    alignment_df["aa_substitutions"] = alignment_df.seq.apply(
#        lambda seq: get_mutations(dms_wt_seq, seq, set(model.mutations))
#    )
#    alignment_df.to_csv(output_df, index=False)
#
#    ret_json = {
#        "generated_by": {"program": "dmsa-pred"},
#        "nodes": defaultdict(dict)
#    }
#    for idx, row in alignment_df.iterrows():
#        ret_json["nodes"][row.strain][f"dms_aa_subs_wrt_{dms_wt_seq_id}"] = row.aa_substitutions
#
#    write_json(ret_json, json_output)


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
@group_options(alignment, dms_wt_seq_id, mut_effects_df, mut_effect_col, experiment_label, output_json, output_df)
def polyclonal_escape_prediction(
    alignment,
    dms_wt_seq_id,
    mut_effects_df,
    mut_effect_col,
    activity_wt_df,
    concentrations,
    experiment_label,
    output_json,
    output_df
):
    """
    Use polyclonal to predict the probability of escape for variants in an
    alignment
    """
    
    # Read in a dataframe with mutational effects. For polyclonal, these
    # are mutational effects in the latent space.
    mut_effects_df = pd.read_csv(mut_effects_df)
    mut_effects_df.rename({mut_effect_col:"escape"}, axis=1, inplace=True)

    # Instantiate a Polyclonal object with the above values and
    # a dataframe of wildtype activities
    model = polyclonal.Polyclonal(
        activity_wt_df=pd.read_csv(activity_wt_df),
        mut_escape_df=mut_effects_df,
        data_to_fit=None,
        sites=tuple(natsort.natsorted(mut_effects_df.site, alg=natsort.ns.SIGNED)),
        alphabet=polyclonal.alphabets.AAS_WITHSTOP_WITHGAP,
    )

    # Read in the input alignment FASTA. Then, for each sequence in the
    # alignment, make a list of all mutations relative to the DMS WT sequence
    alignment_df = fasta_to_df(alignment)
    dms_wt_seq = alignment_df.loc[dms_wt_seq_id, "seq"]
    alignment_df.reset_index(inplace=True)
    #allowed_subs = set(mut_effects_df['aa_substitution'])
    alignment_df["aa_substitutions"] = alignment_df['seq'].apply(
        lambda seq: get_mutations(dms_wt_seq, seq, model.mutations)
    )


    # For each sequence in the alignment, use polyclonal to predict
    # its phenotype given its mutations. Do this for each input
    # concentration.
    concentrations = [float(item) for item in concentrations.split(',')]
    alignment_df = model.prob_escape(
        variants_df=alignment_df,
        concentrations=concentrations
    ).drop("seq", axis=1).reset_index()

    # Write the dataframe of mutations and predicted scores to an
    # output file
    alignment_df['json_label'] = alignment_df['concentration'].apply(
        lambda x: f"{experiment_label}_prob_escape_c_{x}"
    )
    alignment_df.to_csv(output_df, index=False)
    
    # Write the results to an output JSON for use in a Nextstrain workflow
    write_output_json(alignment_df, 'predicted_prob_escape', output_json)

@cli.command(name="escape-fraction")
@option(
    "--condition",
    required=True,
    type=str,
    help="",
)
@group_options(alignment, dms_wt_seq_id, mut_effects_df, mut_effect_col, site_col, experiment_label, output_json, output_df)
def escape_fraction_prediction(
    alignment,
    dms_wt_seq_id,
    mut_effects_df,
    mut_effect_col,
    site_col,
    condition,
    experiment_label,
    output_json,
    output_df
):
    """
    Predict the escape fraction of variants in an alignment
    """

    # Read in data on mutational effects, and subset to data for a specific
    # selection condition
    mut_effects_df = pd.read_csv(mut_effects_df)
    
    # TODO make "condition" an optional argument, and only execute the below
    # line if the argument is provided
    mut_effects_df = mut_effects_df[mut_effects_df['condition'] == condition]

    # Add a column giving the amino-acid substitution for each row
    mut_effects_df['aa_substitution'] = \
        mut_effects_df['wildtype'] + \
        mut_effects_df[site_col].astype('string') + \
        mut_effects_df['mutation'] 

    # Compute the fraction that does NOT escape. These values will be multiplied
    # together below.
    mut_effects_df['non_escape_frac'] = 1 - mut_effects_df[mut_effect_col]
    
    # Read in the input alignment FASTA. Then, for each sequence in the
    # alignment, make a list of all mutations relative to the DMS WT sequence
    alignment_df = fasta_to_df(alignment)
    dms_wt_seq = alignment_df.loc[dms_wt_seq_id, "seq"]
    alignment_df.reset_index(inplace=True)
    allowed_subs = set(mut_effects_df['aa_substitution'])
    alignment_df["aa_substitutions"] = alignment_df['seq'].apply(
        lambda seq: get_mutations(dms_wt_seq, seq, allowed_subs)
    )

    # For each sequence in the alignment, compute its predicted fraction
    # escape based on its mutations
    def predict_escape_fraction(aa_subs):
        data = mut_effects_df[
            mut_effects_df['aa_substitution'].isin(aa_subs.split())
        ]
        return 1 - data['non_escape_frac'].prod()
    
    alignment_df["pred_escape_frac"] = alignment_df['aa_substitutions'].apply(
        lambda x: predict_escape_fraction(x)
    )

    # Write the dataframe of mutations and predicted scores to an
    # output file
    alignment_df['json_label'] = f"{experiment_label}_escape_frac"
    alignment_df.to_csv(output_df, index=False)
    
    # Write the results to an output JSON for use in a Nextstrain workflow
    write_output_json(alignment_df, 'pred_escape_frac', output_json)
    
@cli.command(name="additive-phenotype")
@group_options(alignment, dms_wt_seq_id, mut_effects_df, mut_effect_col, site_col, experiment_label, output_json, output_df)
def additive_phenotype_prediction(
    alignment,
    dms_wt_seq_id,
    mut_effects_df,
    mut_effect_col,
    site_col,
    experiment_label,
    output_json,
    output_df
):
    """
    Predict an additive phenotype of variants in an alignment
    """

    # Read in dataframe with mutation effects
    mut_effects_df = pd.read_csv(mut_effects_df)

    # Add a column that gives the aa substitution
    mut_effects_df['aa_substitution'] = \
        mut_effects_df['wildtype'] + \
        mut_effects_df[site_col].astype('string') + \
        mut_effects_df['mutation']

    # Read in the input alignment FASTA. Then, for each sequence in the
    # alignment, make a list of all mutations relative to the DMS WT sequence
    alignment_df = fasta_to_df(alignment)
    dms_wt_seq = alignment_df.loc[dms_wt_seq_id, "seq"]
    alignment_df.reset_index(inplace=True)
    allowed_subs = set(mut_effects_df['aa_substitution'])
    alignment_df["aa_substitutions"] = alignment_df['seq'].apply(
        lambda seq: get_mutations(dms_wt_seq, seq, allowed_subs)
    )

    # For each sequence in the alignment, compute its predicted phenotype
    # based on its mutations, assuming mutational effects are additive
    def predict_additive_phenotype(aa_subs):
        data = mut_effects_df[
            mut_effects_df['aa_substitution'].isin(aa_subs.split())
        ]
        return data[mut_effect_col].sum()

    alignment_df["pred_score"] = alignment_df["aa_substitutions"].apply(
        lambda x: predict_additive_phenotype(x)
    )

    # Write the dataframe of mutations and predicted scores to an
    # output file
    alignment_df['json_label'] = f"{experiment_label}_pred_score"
    alignment_df.to_csv(output_df, index=False)
    
    # Write the results to an output JSON for use in a Nextstrain workflow
    write_output_json(alignment_df, 'pred_score', output_json)

if __name__ == '__main__':
    cli()
