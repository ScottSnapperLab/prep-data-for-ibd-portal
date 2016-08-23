"""Describe here what this rule accomplishes."""
import logging
logging.basicConfig(filename=snakemake.log.path,level=logging.INFO)

import pandas as pd


# functions
def do_blabla():
    """Do the blabla."""
    pass


# settings
setting1 = snakemake.params.setting1

# input
in_file1 = snakemake.input.in_file1
in_file2 = snakemake.input.in_file2


# output
out_file1 = snakemake.output.out_file1
out_file2 = snakemake.output.out_file2



# business


if setting1:
    # Do something
    pass
else:
    # Do another thing
    pass
