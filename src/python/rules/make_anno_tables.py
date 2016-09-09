"""Extract columns listed in variable ``nmap`` and write out annotation table as CSV."""
import logging
import subprocess as sbp
import shlex as lex

import pandas as pd

from src.python.functions import bool_to_int

# try:
#     snakemake = snakemake
# except NameError:
#     snakemake = None







######## functions

def extract_columns(xls, nmap, converters=None):
    """Return a recoded dataframe with the info needed to add new INFO to the corresponding VCF.

    Columns defined as keys in ``nmap`` will be extracted from each sheet, combined, and renamed.
    """
    assert isinstance(converters, dict) or converters is None

    xls = str(xls)

    sheets = pd.read_excel(xls, sheetname=None, header=0, converters=converters, parse_cols=4)

    # combine into a single table
    df = pd.concat([df for df in sheets.values()], axis=0)

    # change the column names and keep only those in nmap
    df = df.rename(columns=nmap)[list(nmap.values())].copy()

    # break up Chr and Pos and drop Chr_Pos
    df["CHROM"] = df["Chr_Pos"].apply(lambda v: v.split(':')[0])
    df["POS"] = df["Chr_Pos"].apply(lambda v: int(v.split(':')[1]))
    df = df.drop("Chr_Pos", axis=1)

    # Reorder and sort by CHROM and POS
    chrom_pos = ["CHROM","POS"]
    reordered_col = chrom_pos + sorted(list(df.drop(chrom_pos, axis=1).columns.values))
    df = df[reordered_col]
    df = df.sort_values(by=chrom_pos)

    return df


######## constants
convsns = {'BAM CLEAR': bool_to_int,}





def main():
    """Do the main logic."""
    log_form = '%(asctime)s\t%(levelname)s\t%(name)s\t%(message)s'

    # TODO: Fix so that each XLS file's logs are grouped together instead of interleaved
    logging.basicConfig(filename=snakemake.log.path,
                        format=log_form,
                        datefmt='%Y-%m-%dT%H:%M:%S',
                        level=logging.INFO)

    fam_name = snakemake.wildcards.fam_name

    exe = logging.getLogger(fam_name)
    exe.info("-- BEGIN make_anno_tables_1. --")

    ######## params
    # Columms we want mapped to new names we give them
    nmap = snakemake.params.nmap
    # Lines to add to the VCF headers that define the data in the annotation table we will write
    headers = snakemake.params.headers

    ######## input

    xls = snakemake.input.xls


    ######## output

    anno_table_path = snakemake.output.anno_table_path
    anno_headers_path = snakemake.output.anno_headers_path
    # anno_col_name_path = snakemake.output.anno_col_name_path


    ######## business

    exe.info("Converting, extracting and combining columns as needed.")
    df = extract_columns(xls=xls, nmap=nmap, converters=convsns)

    exe.info("Writing table to file.")
    df.to_csv(anno_table_path, sep='\t', index=False)
    #
    # with open(anno_col_name_path,'w') as col_names:
    #     col_names.write(','.join(df.columns.values))

    exe.info("Writing headers to file.")
    with open(anno_headers_path, 'w') as headers_out:
        for h in headers:
            headers_out.write("{header_line}\n".format(header_line=h))


    exe.info("-- END make_anno_tables_1. --")



main()
# if __name__ == '__main__':
#     main()
