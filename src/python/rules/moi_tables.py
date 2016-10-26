#!/usr/bin/env python
"""Produce CSV table to annotate each mode of inherit."""

# Imports
import pandas as pd

from src.python import errors as e

# Metadata
__author__ = "Gus Dunn"
__email__ = "w.gus.dunn@gmail.com"


def load_combine_anno_tables(anno_tables):
    """Return a combined dataframe of all anno_tables."""
    tables = [pd.read_csv(filepath_or_buffer=t, sep='\t') for t in anno_tables]
    return pd.concat(tables)

def validate_mois(annos, mois):
    """Raise error if any encountered MOIs are un-recognized."""
    seen_mois = set(annos.MOI.unique())
    valid_mois = set(mois)

    extra_mois = seen_mois - valid_mois

    if extra_mois:
        msg = "Encountered the following unexpected MOI-values: {extra}. \n\tExpected values are: {valid}"
        raise e.ValidationError(msg.format(extra=list(extra_mois),
                                           valid=list(valid_mois))
                                )

# Functions
def main():
    """Run main work."""
    mois = snakemake.params.mois
    out_dir = snakemake.params.out_dir
    anno_tables = snakemake.input.anno_tables
    moi_tables = snakemake.output.moi_tables

    annos = load_combine_anno_tables(anno_tables=anno_tables)

    validate_mois(annos=annos, mois=mois)

    moi_dfs = {name: pd.DataFrame(df)  for name, df in annos.groupby('MOI')}

    for name, df in moi_dfs.items():
        opath = "{out_dir}/{moi}.csv".format(out_dir=out_dir, moi=name)

        df.to_csv(opath, index=False)



main()
