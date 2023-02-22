
"""
@File: dmsa-pred_cli.py

@Author: Jared Galloway

Command line interface (CLI) for dmsa-pred.
"""

# built-in
import gzip
import glob

# dependencies
import pandas as pd
from click import Path, group, option, argument
import click

# local
# from dmsa-pred import utils
# from dmsa-pred.string import string_ds

import sys
import lzma
from collections import defaultdict
import json
import argparse

# existing deps
import pandas as pd
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

    Here we present a few commands that allow users to
    slice, transform, normalize, fit models to, and more given
    a binary pickle dump'd xarray, usually as a result of running the
    PhIP-Flow pipeline.

    For more information and example workflows please refer to
    the full documentation
    at https://matsengrp.github.io/dmsa-pred/
    """
    pass


@cli.command(name="load-from-csv")
@option(
    "-s",
    "--sample_table",
    required=True,
    type=Path(exists=True),
    help="Path to sample table csv.",
)
@option(
    "-p",
    "--peptide_table",
    required=True,
    type=Path(exists=True),
    help="Path to peptide table csv.",
)
@option(
    "-c",
    "--counts_matrix",
    required=True,
    type=Path(exists=True),
    help="Path to counts matrix csv.",
)
@option(
    "-o",
    "--output",
    required=True,
    help="Path where the phip dataset will be dump'd to netCDF",
)
def load_from_csv(
    sample_table,
    peptide_table,
    counts_matrix,
    output,
):
    """
    Load and dump xarray dataset given a set of wide csv's

    Using this command usually means you have either:

    1. Decided to store the output of your analysis in the form
       of wide csv's instead of a pickle dump'd binary for
       longer-term storage.

    2. Created your own enrichment matrix
       without the help of the phip-flow alignment pipeline.

    \f

    .. note::
      In the case of #2, please note that your matrix data
      must be numeric and have shape (len(peptide_table), len(sample_table)).
      Finally, you must include
      :ref:`pipeline outputs <sec_pipeline_outputs>`

    .. note::
      Currently only accepting a single enrichment matrix.
    """

    ds = utils.dataset_from_csv(counts_matrix, peptide_table, sample_table)
