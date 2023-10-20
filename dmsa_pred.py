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
import numpy as np
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

def get_mutations(seq1, seq2):
    """Make a list of amino-acid mutations between two sequences"""

    assert len(seq1) == len(seq2)
    return [
        f"{aa1}{pos+1}{aa2}"
        for pos, (aa1, aa2) in enumerate(zip(seq1, seq2))
        if aa1 != aa2
    ]

def parse_aa_substitutions_from_alignment(
        alignment,
        dms_wt_seq_id,
        muts_with_measurements,
        allow_aa_subs_at_unmeasured_sites,
        allow_unmeasured_aa_subs_at_these_sites
    ):
    """
    Parse amino-acid substitutions from an input alignment
    """

    # Read in the input alignment FASTA. Then, for each sequence in the
    # alignment, make a list of all mutations relative to the DMS WT sequence
    alignment_df = fasta_to_df(alignment)
    dms_wt_seq = alignment_df.loc[dms_wt_seq_id, "seq"]
    alignment_df.reset_index(inplace=True)
    alignment_df["all_aa_substitutions"] = alignment_df['seq'].apply(
        lambda seq: get_mutations(dms_wt_seq, seq)
    )

    # If indicated, make a list of sites from the alignment that completely lack DMS data (i.e.
    # zero mutation effects were measured at those sites), and use it below to avoid masking
    # sequences with substitutions at those sites
    unmeasured_sites_to_allow = []
    if allow_aa_subs_at_unmeasured_sites:
        sites_with_measurements = set([mut[1:-1] for mut in muts_with_measurements])
        unmeasured_sites_to_allow += [
            str(site) for (site, aa) in enumerate(dms_wt_seq, 1)
            if str(site) not in sites_with_measurements
        ]

    # If provided, read in a user-provided list of sites to add to unmeasured_sites_to_allow
    if allow_unmeasured_aa_subs_at_these_sites != 'None':
        with open(allow_unmeasured_aa_subs_at_these_sites) as f:
            unmeasured_sites_to_allow += f.read().split('\n')

    # For each sequence, make a list of mutations with measurements, a list of
    # unmeasured mutations, and a list of unmeasured mutations that are disallowed
    alignment_df["measured_aa_substitutions"] = alignment_df["all_aa_substitutions"].apply(
        lambda aa_substitutions: [
            sub for sub in aa_substitutions
            if sub in muts_with_measurements
        ]
    )
    alignment_df["unmeasured_aa_substitutions"] = alignment_df["all_aa_substitutions"].apply(
        lambda aa_substitutions: [
            sub for sub in aa_substitutions
            if sub not in muts_with_measurements
        ]
    )
    alignment_df["disallowed_aa_substitutions"] = alignment_df["unmeasured_aa_substitutions"].apply(
        lambda aa_substitutions: [
            sub for sub in aa_substitutions
            if sub[1:-1] not in unmeasured_sites_to_allow
        ]
    )
    alignment_df["n_disallowed_aa_substitutions"] = alignment_df["disallowed_aa_substitutions"].apply(
        lambda x: len(x)
    )

    # For each column listing aa substitutions, turn lists into space-delimited strings. Also,
    # remove the sequence column to save space.
    aa_subs_cols = [
        'all_aa_substitutions', 'measured_aa_substitutions', 'unmeasured_aa_substitutions',
        'disallowed_aa_substitutions'
        ]
    for col in aa_subs_cols:
        alignment_df[col] = alignment_df[col].apply(lambda x: ' '.join(x))
    del alignment_df['seq']

    return alignment_df

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
mutation_col = option(
    "--mutation-col",
    required=True,
    type=str,
    help="a column from mut_effects_df that gives the genotype of a variant i.e. <wt><site><mutation>",
)
mask_seqs_with_disallowed_aa_subs = option(
    "--mask-seqs-with-disallowed-aa-subs",
    required=False,
    default=True,
    type=bool,
    help="if True (default), sequences with unmeasured aa substitutions will be excluded from the output JSON. Note: certain aa substitutions are allowed to be missing if allow-aa-subs-at-unmeasured-sites is True or allow-unmeasured-aa-subs-at-these-sites is True. See the documentation of those arguments for more detail.",
)
allow_aa_subs_at_unmeasured_sites = option(
    "--allow-aa-subs-at-unmeasured-sites",
    required=False,
    default=False,
    type=bool,
    help="some sites in the alignment might completely lack DMS data for any mutation at that site. If this argument is set to True, then these sites will not be considered when identifying sequences with unmeasured aa substitutions for masking when mask_seqs_with_disallowed_mutations is True",
)
allow_unmeasured_aa_subs_at_these_sites = option(
    "--allow-unmeasured-aa-subs-at-these-sites",
    required=False,
    default='None',
    type=str,
    help="the path to a file with a list of newline-deliminted sites. If a file path is provided, then the listed sites will not be considered when identifying sequences with unmeasured aa substitutions for masking when mask_seqs_with_disallowed_mutations is True. If allow-aa-subs-at-unmeasured-sites is also True, then the sites from each list are simply combined into one list. A string value of 'None' will be ignored.",
)
min_pred_pheno = option(
    "--min-pred-pheno",
    required=False,
    type=float,
    help="if provided, predicted phenotypes will be clipped at this lower value",
)
max_pred_pheno = option(
    "--max-pred-pheno",
    required=False,
    type=float,
    help="if provided, predicted phenotypes will be clipped at this upper value",
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

# TODO: do we still need this?
@cli.command(name="call-mutations")
@group_options(alignment, dms_wt_seq_id, output_json, output_df)
def call_mutations(
    alignment,
    dms_wt_seq_id,
    output_json,
    output_df
):
    """
    Read in the input alignment FASTA. Then, for each sequence in the
    # alignment, make a list of all mutations relative to the DMS WT sequence
    """

    alignment_df = fasta_to_df(alignment)
    dms_wt_seq = alignment_df.loc[dms_wt_seq_id, "seq"]
    alignment_df.reset_index(inplace=True)
    alignment_df["aa_substitutions"] = alignment_df.seq.apply(
        lambda seq: get_mutations(dms_wt_seq, seq, set(model.mutations))
    )
    alignment_df.to_csv(output_df, index=False)

    ret_json = {
        "generated_by": {"program": "dmsa-pred"},
        "nodes": defaultdict(dict)
    }
    for idx, row in alignment_df.iterrows():
        ret_json["nodes"][row.strain][f"dms_aa_subs_wrt_{dms_wt_seq_id}"] = row.aa_substitutions

    write_json(ret_json, json_output)


@cli.command(name="polyclonal-escape")
@option(
    "--activity-wt-df",
    required=True,
    type=Path(exists=True),
    help="",
)
@option(
    "--concentrations",
    required=False,
    default="0.0",
    type=str,
    help="",
)
@option(
    "--icxx",
    required=False,
    default=0.0,
    type=float,
    help="",
)
@group_options(
    alignment, 
    dms_wt_seq_id, 
    mut_effects_df, 
    mut_effect_col, 
    mutation_col,
    mask_seqs_with_disallowed_aa_subs,
    allow_aa_subs_at_unmeasured_sites,
    allow_unmeasured_aa_subs_at_these_sites,
    min_pred_pheno,
    max_pred_pheno,
    experiment_label, 
    output_json, 
    output_df
)
def polyclonal_escape_prediction(
    alignment,
    dms_wt_seq_id,
    mut_effects_df,
    mut_effect_col,
    mutation_col,
    mask_seqs_with_disallowed_aa_subs,
    allow_aa_subs_at_unmeasured_sites,
    allow_unmeasured_aa_subs_at_these_sites,
    min_pred_pheno,
    max_pred_pheno,
    activity_wt_df,
    concentrations,
    icxx,
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
    mut_effects_df.rename(
        {
            mut_effect_col:"escape",
            mutation_col:"mutation"
        }, 
        axis=1, 
        inplace=True
    )
    muts_with_measurements = set(mut_effects_df['mutation'])

    # Instantiate a Polyclonal object with the above values and
    # a dataframe of wildtype activities
    model = polyclonal.Polyclonal(
        activity_wt_df=pd.read_csv(activity_wt_df),
        mut_escape_df=mut_effects_df,
        data_to_fit=None,
        sites=tuple(natsort.natsorted(mut_effects_df.site, alg=natsort.ns.SIGNED)),
        alphabet=polyclonal.alphabets.AAS_WITHSTOP_WITHGAP,
    )

    # Parse aa substitutions from the alignment
    alignment_df = parse_aa_substitutions_from_alignment(
        alignment,
        dms_wt_seq_id,
        muts_with_measurements,
        allow_aa_subs_at_unmeasured_sites,
        allow_unmeasured_aa_subs_at_these_sites
    )
    alignment_df["aa_substitutions"] = alignment_df['measured_aa_substitutions']

    # Use polyclonal to compute phenotypes of each sequence in the alignment
    node_attrs = []
    if icxx != 0.0:

        # Check that the `concentrations` option is not also specified
        if concentrations != "0.0":
            raise ValueError("the options `icxx` and `concentrations` cannot both be used at the same time")

        # Compute the delta icxx
        col = f"{experiment_label}_IC{int(icxx*100)}"
        model_preds = model.icXX(
            variants_df=alignment_df,
            x=icxx,
            col=col 
        )
        alignment_df = alignment_df.merge(model_preds[["strain",col]], on="strain")
        wt_icxx = model_preds.loc[alignment_df.strain == dms_wt_seq_id, col].values[0]
        alignment_df[f"{col}_log_fold_change"] = np.log10(alignment_df[col] / wt_icxx)
        node_attrs.append(f"{col}_log_fold_change")

        # If indicated, clip scores at an upper/lower bound, while adding a column that
        # gives the unaltered phenotype
        if min_pred_pheno or max_pred_pheno:
            alignment_df[f"raw_{col}_log_fold_change"] = alignment_df[f"{col}_log_fold_change"]
        if isinstance(min_pred_pheno, float):
            alignment_df[f"{col}_log_fold_change"].clip(lower=min_pred_pheno, inplace=True)
        if isinstance(max_pred_pheno, float):
            alignment_df[f"{col}_log_fold_change"].clip(upper=max_pred_pheno, inplace=True)

    elif concentrations != "0.0":

        concentrations = [float(item) for item in concentrations.split(',')]
            
        model_preds = model.prob_escape(
            variants_df=alignment_df,
            concentrations=concentrations
        )

        for con, con_model_preds in model_preds.groupby('concentration'):
            alignment_df = alignment_df.merge(
                con_model_preds[["strain","predicted_prob_escape"]], on="strain"
            )
            col = f"{experiment_label}_prob_escape_c_{con}"
            alignment_df.rename({"predicted_prob_escape":col}, axis=1, inplace=True)
            node_attrs.append(col)

            # If indicated, clip scores at an upper/lower bound, while adding a column that
            # gives the unaltered phenotype
            if min_pred_pheno or max_pred_pheno:
                alignment_df[f"raw_{col}"] = alignment_df[col]
            if isinstance(min_pred_pheno, float):
                alignment_df[col].clip(lower=min_pred_pheno, inplace=True)
            if isinstance(max_pred_pheno, float):
                alignment_df[col].clip(upper=max_pred_pheno, inplace=True)
    else:
        raise ValueError("either the `icxx` option or the `concentrations` option must be provided")

    # Write the dataframe of mutations and predicted scores to an output file
    alignment_df.to_csv(output_df, index=False)

    # If indicated, mask sequences with disallowed aa substitutions by dropping them from
    # the dataframe before writing the JSON
    if mask_seqs_with_disallowed_aa_subs:
        alignment_df = alignment_df[alignment_df['n_disallowed_aa_substitutions'] == 0]

    # Write data to output JSON file
    ret_json = {
        "generated_by": {"program": "dmsa-pred"},
        "nodes": defaultdict(dict)
    }
    for idx, row in alignment_df.iterrows():
        for attr in node_attrs:
            ret_json["nodes"][row.strain][attr] = row[attr]
    write_json(ret_json, output_json)

@cli.command(name="phenotype-prediction")
@option(
    "--model-type",
    required=False,
    default="additive",
    type=str,
    help="",
)
@group_options(
    alignment, 
    dms_wt_seq_id, 
    mut_effects_df, 
    mut_effect_col, 
    mutation_col,
    mask_seqs_with_disallowed_aa_subs,
    allow_aa_subs_at_unmeasured_sites,
    allow_unmeasured_aa_subs_at_these_sites,
    min_pred_pheno,
    max_pred_pheno,
    experiment_label, 
    output_json, 
    output_df
)
def predict_phenotypes(
    model_type,
    alignment,
    dms_wt_seq_id,
    mut_effects_df,
    mut_effect_col,
    mutation_col,
    mask_seqs_with_disallowed_aa_subs,
    allow_aa_subs_at_unmeasured_sites,
    allow_unmeasured_aa_subs_at_these_sites,
    min_pred_pheno,
    max_pred_pheno,
    experiment_label,
    output_json,
    output_df
):
    """
    Predict the phenotype of variants in an alignment using either the additive effects model or a
    model of phenotype by the product of (1-phenotype) fractions.
    """

    # Read in dataframe with mutational effects and get the set of mutations
    # with measurements
    mut_effects_df = pd.read_csv(mut_effects_df)
    muts_with_measurements = set(mut_effects_df[mutation_col])


    # Parse aa substitutions from the alignment
    alignment_df = parse_aa_substitutions_from_alignment(
        alignment,
        dms_wt_seq_id,
        muts_with_measurements,
        allow_aa_subs_at_unmeasured_sites,
        allow_unmeasured_aa_subs_at_these_sites
    )

    # We are currently not using the phenotype fraction, but we'll leave it here
    # for now in case we want to easily impliment it again.

    if model_type == "fraction_escape":
        
        # For each sequence in the alignment, compute its predicted fraction
        # phenotype based on its mutations
        mut_effects_df['non_phenotype_frac'] = 1 - mut_effects_df[mut_effect_col]

        def predict_phenotype_fraction(aa_subs):
            data = mut_effects_df[
                mut_effects_df[mutation_col].isin(aa_subs.split())
            ]
            return 1 - data['non_phenotype_frac'].prod()

        alignment_df['pred_phenotype'] = alignment_df['measured_aa_substitutions'].apply(
            lambda x: predict_phenotype_fraction(x)
        )

    elif model_type == "additive":

        # For each sequence in the alignment, compute its predicted phenotype
        # based on its mutations, assuming mutational effects are additive
        def predict_additive_phenotype(aa_subs):
            data = mut_effects_df[
                mut_effects_df[mutation_col].isin(aa_subs.split())
            ]
            return data[mut_effect_col].sum()

        alignment_df['pred_phenotype'] = alignment_df["measured_aa_substitutions"].apply(
            lambda x: predict_additive_phenotype(x)
        )

    else:
        raise ValueError(f"model_type {model_type} is unknown, please specify either 'additive' or 'fraction_escape'")

    # If indicated, clip scores at an upper/lower bound, while adding a column that
    # gives the unaltered phenotype
    if min_pred_pheno or max_pred_pheno:
        alignment_df['raw_pred_phenotype'] = alignment_df['pred_phenotype']
    if isinstance(min_pred_pheno, float):
        alignment_df['pred_phenotype'].clip(lower=min_pred_pheno, inplace=True)
    if isinstance(max_pred_pheno, float):
        alignment_df['pred_phenotype'].clip(upper=max_pred_pheno, inplace=True)

    # Write the dataframe of mutations and predicted scores to an output file
    # alignment_df['json_label'] = f"{experiment_label}_{model_type}_phenotype"
    alignment_df.to_csv(output_df, index=False)
    
    # If indicated, mask sequences with disallowed aa substitutions by dropping them from
    # the dataframe before writing the JSON
    if mask_seqs_with_disallowed_aa_subs:
        alignment_df = alignment_df[alignment_df['n_disallowed_aa_substitutions'] == 0]

    # Write data to output JSON file
    ret_json = {
        "generated_by": {"program": "dmsa-pred"},
        "nodes": defaultdict(dict)
    }
    for idx, row in alignment_df.iterrows():
        ret_json["nodes"][row.strain][experiment_label] = row['pred_phenotype']
    write_json(ret_json, output_json)


if __name__ == '__main__':
    cli()
