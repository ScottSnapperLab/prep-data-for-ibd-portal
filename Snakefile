"""Snakemake file."""
import os

from pathlib import Path

import yaml

import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

import munch

from src.python.functions import *


def pathify_by_key_ends(dictionary):
    """Return a dict that has had all values with keys marked as '*_PATH' or '*_DIR' converted to Path() instances."""
    for key, value in dictionary.items():
        if isinstance(value, dict):
            pathify_by_key_ends(value)
        elif key.endswith("_PATH") or key.endswith("_DIR"):
            dictionary[key] = Path(value)

    return dictionary


class MyRun(object):

    """Initialize and manage information common to the whole run."""

    def __init__(self, cfg):
        """Initialize common information for a run."""
        assert isinstance(cfg, dict)

        common = cfg["COMMON"]

        self.globals = munch.Munch()
        self.cfg = cfg
        self.name = common["RUN_NAME"]
        self.d = common["SHARED"]
        self.out_dir = Path("{base_dir}/{run_name}".format(base_dir=common["OUT_DIR"],
                                                           run_name=self.name
                                                           )
                            )
        self.log_dir = self.out_dir / "logs"

class MyRule(object):

    """Manage the initialization and deployment of rule-specific information."""

    def __init__(self, run, name):
        """Initialize logs, inputs, outputs, params, etc for a single rule."""
        assert isinstance(run, MyRun)

        self.run = run
        self.name = name.lower()
        self.log_dir = run.log_dir / self.name
        self.log = self.log_dir / "{name}.log".format(name=self.name)
        self.out_dir = run.out_dir / self.name
        self.i = munch.Munch() # inputs
        self.o = munch.Munch() # outputs
        self.p = munch.Munch() # params

        self._import_config_dict()

    def _import_config_dict(self):
        """Inport configuration values set for this rule so they are directly accessable as attributes."""
        try:
            for key, val in self.run.cfg[self.name.upper()].items():
                self.__setattr__(key, val)
            self.cfg = True
        except KeyError:
            self.cfg = False



#### COMMON RUN STUFF ####
ORIGINAL_CONFIG_AS_STRING = yaml.dump(config, default_flow_style=False)
config = pathify_by_key_ends(config)
config = munch.munchify(config)

RUN = MyRun(cfg=config)

PRE = []
VALIDATE_INPUT = []
CONVERSIONS = []
MERGES = []

# add specific useful stuff to RUN
RUN.globals.fam_names = [vcf.stem for vcf in config.VALIDATE_INPUT_VCFS.IN.VCF_DIR.glob("*.vcf")]


############ BEGIN PIPELINE RULES ############
# ------------------------- #
#### SAVE_RUN_CONFIG ####
SAVE_RUN_CONFIG = MyRule(run=RUN, name="SAVE_RUN_CONFIG")
SAVE_RUN_CONFIG.o.file = RUN.out_dir / "{NAME}.yaml".format(NAME=RUN.name)



rule save_run_config:
    priority: 50
    input:
    output:
        file=str(SAVE_RUN_CONFIG.o.file)

    run:
        with open(output.file, 'w') as cnf_out:
            cnf_out.write(ORIGINAL_CONFIG_AS_STRING)

PRE.append(rules.save_run_config.output)
CONVERSIONS.append(rules.save_run_config.output)
MERGES.append(rules.save_run_config.output)












# ------------------------- #
#### CREATE_SEQ_DICT ####
CREATE_SEQ_DICT = MyRule(run=RUN, name="CREATE_SEQ_DICT")

CREATE_SEQ_DICT.o.fasta_dict = RUN.d["REFERENCE_FASTA_PATH"].parent / "{stem}.dict".format(stem=RUN.d["REFERENCE_FASTA_PATH"].stem)

# ---
rule create_seq_dict:
    priority: 99
    log:
        path=str(CREATE_SEQ_DICT.log)

    input:
        ref_fasta=str(RUN.d["REFERENCE_FASTA_PATH"]),

    output:
        fasta_dict=str(CREATE_SEQ_DICT.o.fasta_dict),

    shell:
        "picard CreateSequenceDictionary "
        "R={input.ref_fasta} "
        "O={output.fasta_dict} "
        "&> {log.path} "

PRE.append(rules.create_seq_dict.output)



# ------------------------- #
#### CREATE_SEQ_FAI ####
CREATE_SEQ_FAI = MyRule(run=RUN, name="CREATE_SEQ_FAI")

CREATE_SEQ_FAI.o.fasta_fai = "{fas}.fai".format(fas=RUN.d["REFERENCE_FASTA_PATH"])


# ---
rule create_seq_fai:
    priority: 98
    log:
        path=str(CREATE_SEQ_FAI.log)

    input:
        ref_fasta=str(RUN.d["REFERENCE_FASTA_PATH"]),

    output:
        fasta_fai=str(CREATE_SEQ_FAI.o.fasta_fai),

    shell:
        "samtools faidx "
        "{input.ref_fasta} "
        "&> {log.path} "

PRE.append(rules.create_seq_fai.output)




# ------------------------- #
#### VALIDATE_INPUT_VCFS ####
VALIDATE_INPUT_VCFS = MyRule(run=RUN, name="VALIDATE_INPUT_VCFS")

# input
VALIDATE_INPUT_VCFS.i.input_vcf = str(VALIDATE_INPUT_VCFS.IN.VCF_DIR / "{fam_name}.vcf")

# output
VALIDATE_INPUT_VCFS.o.sentinel_wldcd = str(VALIDATE_INPUT_VCFS.out_dir / "{fam_name}_validated_on.txt" )

VALIDATE_INPUT_VCFS.o.sentinels_expanded = expand(VALIDATE_INPUT_VCFS.o.sentinel_wldcd, fam_name=RUN.globals.fam_names)


# ---
rule validate_input_vcfs:
    priority: 90
    log:
        path=str(VALIDATE_INPUT_VCFS.log)

    input:
        rules.create_seq_dict.output,
        rules.create_seq_fai.output,
        input_vcf=VALIDATE_INPUT_VCFS.i.input_vcf,
        ref_fasta=str(RUN.d["REFERENCE_FASTA_PATH"]),


    output:
        sentinels=VALIDATE_INPUT_VCFS.o.sentinel_wldcd,

    shell:
        "gatk -T ValidateVariants "
        "-R {input.ref_fasta} "
        "-V {input.input_vcf} "
        "--warnOnErrors "
        # "--validationTypeToExclude ALL "
        "&> {log.path} "
        "&& echo $(date) > {output.sentinels}"

CONVERSIONS.append(VALIDATE_INPUT_VCFS.o.sentinels_expanded)


# ------------------------- #
#### RECODE_INPUT_VCFS ####
RECODE_INPUT_VCFS = MyRule(run=RUN, name="RECODE_INPUT_VCFS")

# params

# input
RECODE_INPUT_VCFS.i.vcf = str(VALIDATE_INPUT_VCFS.IN.VCF_DIR / "{fam_name}.vcf")

# output
RECODE_INPUT_VCFS.o.vcf_wldcd = str(RECODE_INPUT_VCFS.out_dir / "{fam_name}.vcf")

RECODE_INPUT_VCFS.o.vcf_expanded = expand(RECODE_INPUT_VCFS.o.vcf_wldcd, fam_name=RUN.globals.fam_names)

# ---
rule recode_input_vcfs:
    log:
        path=str(RECODE_INPUT_VCFS.log_dir / "{fam_name}.log")

    params:
        # param_1=param_1,

    input:
        validation_complete=VALIDATE_INPUT_VCFS.o.sentinel_wldcd,
        vcf=RECODE_INPUT_VCFS.i.vcf,

    output:
        vcf=RECODE_INPUT_VCFS.o.vcf_wldcd,

    shell:
        "cat {input.vcf} | "
        "sed -u -e 's/AD,Number=./AD,Number=R/g' "
        "-e 's/MLEAC,Number=A/MLEAC,Number=R/g' "
        "-e 's/MLEAF,Number=A/MLEAF,Number=R/g' "
        "-e 's/PL,Number=./PL,Number=G/g' > "
        "{output.vcf} 2> {log.path}"

CONVERSIONS.append(RECODE_INPUT_VCFS.o.vcf_expanded)


# ------------------------- #
#### BGZIP_AND_IDX_VCFS ####
BGZIP_AND_IDX_VCFS = MyRule(run=RUN, name="BGZIP_AND_IDX_VCFS")

# params

# input
BGZIP_AND_IDX_VCFS.i.vcf_wldcd = RECODE_INPUT_VCFS.o.vcf_wldcd

# output
BGZIP_AND_IDX_VCFS.o.vcf_bgz_wldcd = str(BGZIP_AND_IDX_VCFS.out_dir / "{fam_name}.vcf.bgz")
BGZIP_AND_IDX_VCFS.o.vcf_bgz_tbi_wldcd = str(BGZIP_AND_IDX_VCFS.out_dir / "{fam_name}.vcf.bgz.tbi")

BGZIP_AND_IDX_VCFS.o.vcf_bgz_expanded = expand(BGZIP_AND_IDX_VCFS.o.vcf_bgz_wldcd, fam_name=RUN.globals.fam_names)
BGZIP_AND_IDX_VCFS.o.vcf_bgz_tbi_expanded = expand(BGZIP_AND_IDX_VCFS.o.vcf_bgz_tbi_wldcd, fam_name=RUN.globals.fam_names)

# ---
rule bgzip_and_idx_vcfs:
    log:
        path=str(BGZIP_AND_IDX_VCFS.log_dir / "{fam_name}")

    params:

    input:
        vcf=BGZIP_AND_IDX_VCFS.i.vcf_wldcd,

    output:
        vcf_bgz=BGZIP_AND_IDX_VCFS.o.vcf_bgz_wldcd,
        vcf_bgz_tbi=BGZIP_AND_IDX_VCFS.o.vcf_bgz_tbi_wldcd,

    shell:
        "bgzip -c {input.vcf} > {output.vcf_bgz} 2>> {log.path}.bgzip.log &&"
        "tabix -p vcf {output.vcf_bgz} > {log.path}.tabix.log 2>&1"

CONVERSIONS.append(BGZIP_AND_IDX_VCFS.o.vcf_bgz_expanded + \
                 BGZIP_AND_IDX_VCFS.o.vcf_bgz_tbi_expanded)






# ------------------------- #
#### MAKE_ANNO_TABLES ####
MAKE_ANNO_TABLES = MyRule(run=RUN, name="MAKE_ANNO_TABLES")

# params
# :: see rule definition below

# input
# :: see rule definition below

# output
MAKE_ANNO_TABLES.o.anno_table_path_wldcd = str(MAKE_ANNO_TABLES.out_dir / "{fam_name}.csv")
MAKE_ANNO_TABLES.o.anno_table_path_bgz_wldcd = str(MAKE_ANNO_TABLES.out_dir / "{fam_name}.csv.gz")
MAKE_ANNO_TABLES.o.anno_table_path_tbi_wldcd = str(MAKE_ANNO_TABLES.out_dir / "{fam_name}.csv.gz.tbi")
MAKE_ANNO_TABLES.o.anno_headers_path_wldcd = str(MAKE_ANNO_TABLES.out_dir / "{fam_name}.hdr")
# MAKE_ANNO_TABLES.o.anno_col_name_path_wldcd = str(MAKE_ANNO_TABLES.out_dir / "{fam_name}.cols")

# :: determine the names of output files so snakemake can do its dependency checks
# XLS_NAMES = list(MAKE_ANNO_TABLES.IN.XLS_DIR.glob("*.xlsx"))
MAKE_ANNO_TABLES.o.anno_table_path_expanded = expand(MAKE_ANNO_TABLES.o.anno_table_path_wldcd, fam_name=RUN.globals.fam_names)
MAKE_ANNO_TABLES.o.anno_table_path_bgz_expanded = expand(MAKE_ANNO_TABLES.o.anno_table_path_bgz_wldcd, fam_name=RUN.globals.fam_names)
MAKE_ANNO_TABLES.o.anno_table_path_tbi_expanded = expand(MAKE_ANNO_TABLES.o.anno_table_path_bgz_wldcd, fam_name=RUN.globals.fam_names)
MAKE_ANNO_TABLES.o.anno_headers_path_expanded = expand(MAKE_ANNO_TABLES.o.anno_headers_path_wldcd, fam_name=RUN.globals.fam_names)
# MAKE_ANNO_TABLES.o.anno_col_name_path_expanded = expand(MAKE_ANNO_TABLES.o.anno_col_name_path_wldcd, fam_name=RUN.globals.fam_names)

MAKE_ANNO_TABLES.OUTPUT_1 = MAKE_ANNO_TABLES.o.anno_table_path_expanded + \
                            MAKE_ANNO_TABLES.o.anno_headers_path_expanded #+\
                            # MAKE_ANNO_TABLES.o.anno_col_name_path_expanded


MAKE_ANNO_TABLES.OUTPUT_2 =   MAKE_ANNO_TABLES.o.anno_table_path_bgz_expanded + \
                              MAKE_ANNO_TABLES.o.anno_table_path_tbi_expanded

# ---
rule make_anno_tables_1:
    log:
        path=str(MAKE_ANNO_TABLES.log_dir / "{fam_name}.log")

    params:
        nmap=MAKE_ANNO_TABLES.PARAMS.NMAP,
        headers=MAKE_ANNO_TABLES.PARAMS.HEADERS,

    input:
        xls=str(MAKE_ANNO_TABLES.IN.XLS_DIR / "{fam_name}.xls"),

    output:
        anno_table_path=MAKE_ANNO_TABLES.o.anno_table_path_wldcd,
        anno_headers_path=MAKE_ANNO_TABLES.o.anno_headers_path_wldcd,
        # anno_col_name_path=MAKE_ANNO_TABLES.o.anno_col_name_path_wldcd



    script:
        "src/python/rules/make_anno_tables.py"

CONVERSIONS.append(MAKE_ANNO_TABLES.OUTPUT_1)

# ---
rule make_anno_tables_2:
    log:
        path=str(MAKE_ANNO_TABLES.log_dir / "{fam_name}")
    #
    # params:
    #     nmap=MAKE_ANNO_TABLES.PARAMS.NMAP,
    #     headers=MAKE_ANNO_TABLES.PARAMS.HEADERS,
    input:
        anno_table=MAKE_ANNO_TABLES.o.anno_table_path_expanded,

    output:
        anno_table_path_bgz=MAKE_ANNO_TABLES.o.anno_table_path_bgz_wldcd,
        anno_table_path_tbi=MAKE_ANNO_TABLES.o.anno_table_path_tbi_wldcd,


    shell:
        "bgzip -c {input.anno_table} > {output.anno_table_path_bgz} 2>> {log.path}.bgzip.log && "
        "tabix -s1 -b2 -e2 -S1 {output.anno_table_path_bgz} > {log.path}.tabix.log 2>&1"

CONVERSIONS.append(MAKE_ANNO_TABLES.OUTPUT_2)


# ------------------------- #
#### ANNOTATE_VCFS ####
ANNOTATE_VCFS = MyRule(run=RUN, name="ANNOTATE_VCFS")

# params
ANNOTATE_VCFS.p.cols = ANNOTATE_VCFS.PARAMS.COLS

# input
ANNOTATE_VCFS.i.vcf_wldcd = str(BGZIP_AND_IDX_VCFS.out_dir / '{fam_name}.vcf.bgz')
# ANNOTATE_VCFS.i.vcf_wldcd = str(ANNOTATE_VCFS.IN.VCF_DIR / '{fam_name}.vcf')
ANNOTATE_VCFS.i.anno_table_path_bgz_wldcd = MAKE_ANNO_TABLES.o.anno_table_path_bgz_wldcd
ANNOTATE_VCFS.i.anno_headers_path_wldcd = MAKE_ANNO_TABLES.o.anno_headers_path_wldcd


# output
ANNOTATE_VCFS.o.bcf_wldcd = str(ANNOTATE_VCFS.out_dir / '{fam_name}.bcf')

ANNOTATE_VCFS.o.bcf_expanded = expand(ANNOTATE_VCFS.o.bcf_wldcd, fam_name=RUN.globals.fam_names)

# ---
rule annotate_vcfs:
    log:
        path=str(ANNOTATE_VCFS.log_dir / "{fam_name}.log")

    params:
        cols=ANNOTATE_VCFS.p.cols,

    input:
        vcf=ANNOTATE_VCFS.i.vcf_wldcd,
        anno_table_path_bgz=ANNOTATE_VCFS.i.anno_table_path_bgz_wldcd,
        anno_headers_path=ANNOTATE_VCFS.i.anno_headers_path_wldcd,

    output:
        bcf=ANNOTATE_VCFS.o.bcf_wldcd,

    shell:
        "bcftools annotate "
        "-x INFO "
        # "-a {input.anno_table_path_bgz} "
        # "-h {input.anno_headers_path} "
        # "-c {params.cols} "
        "{input.vcf} "
        "-Ou "
        "-o {output.bcf} > {log.path} 2>&1"

CONVERSIONS.append(ANNOTATE_VCFS.o.bcf_expanded)


# ------------------------- #
#### LEFT_ALN_AND_NORM ####
LEFT_ALN_AND_NORM = MyRule(run=RUN, name="LEFT_ALN_AND_NORM")

# params
LEFT_ALN_AND_NORM.p.ref_fasta = str(RUN.d["REFERENCE_FASTA_PATH"])

# input
LEFT_ALN_AND_NORM.i.bcf_wldcd = str(ANNOTATE_VCFS.out_dir / "{fam_name}.bcf")

# output
LEFT_ALN_AND_NORM.o.vcf_wldcd = str(LEFT_ALN_AND_NORM.out_dir / "{fam_name}.vcf")

LEFT_ALN_AND_NORM.o.vcf_expanded = expand(LEFT_ALN_AND_NORM.o.vcf_wldcd, fam_name=RUN.globals.fam_names)


# ---
rule left_aln_and_norm:
    log:
        path=str(LEFT_ALN_AND_NORM.log_dir / "{fam_name}")

    params:
        ref_fasta=LEFT_ALN_AND_NORM.p.ref_fasta,

    input:
        bcf =LEFT_ALN_AND_NORM.i.bcf_wldcd ,

    output:
        vcf=LEFT_ALN_AND_NORM.o.vcf_wldcd,

    shell:
        "vt normalize {input.bcf} "
        "-r {params.ref_fasta} "
        "-o + 2> {log.path}.vt.normalize.log  | "
        "vt decompose - "
        "-s "
        "-o {output.vcf} > {log.path}.vt.decompse.log 2>&1"



CONVERSIONS.append(LEFT_ALN_AND_NORM.o.vcf_expanded)


# # ------------------------- #
# #### MERGE_BCFS ####
# MERGE_BCFS = MyRule(run=RUN, name="MERGE_BCFS")
#
# # params
# MERGE_BCFS.p.param_1 = MERGE_BCFS.PARAMS.param_1
#
# # input
# MERGE_BCFS.i.input_1 = str(MERGE_BCFS.IN.input_1 / "{something}.ext")
#
# # output
# MERGE_BCFS.o.output_1 = str(MERGE_BCFS.out_dir / "{something}.ext")
#
# # ---
# rule merge_bcfs:
#     log:
#         path=str(MERGE_BCFS.log)
#
#     params:
#         param_1=MERGE_BCFS.p.param_1,
#
#     input:
#         input_1=MERGE_BCFS.i.input_1,
#
#     output:
#         output_1=MERGE_BCFS.o.output_1,
#
#     shell:
#         ""
#
# MERGES.append(rules.merge_bcfs.output)





# ------------------------- #
#### PRE ####
# ---
rule pre:
    input: PRE


# ------------------------- #
#### CONVERSIONS ####
# ---
rule conversions:
    input: CONVERSIONS


# ------------------------- #
#### MERGES ####
# ---
rule merges:
    input: MERGES


# ------------------------- #
#### CLEAN_CONVERSIONS ####
# ---
rule clean_conversions:
    shell:
        "rm -rf {RUN.out_dir} && find data/raw -type f -name '*.idx' -print0 | xargs -0 rm &> /dev/null || true"

# ------------------------- #
#### CLEAN_PRE ####
# ---
rule clean_pre:
    shell:
        "find data/external -type f -name '*.dict' -print0 | xargs -0 rm &> /dev/null || true && "
        "find data/external -type f -name '*.fai' -print0 | xargs -0 rm &> /dev/null || true"

# ------------------------- #
#### CLEAN_ALL ####
# ---
rule clean_all:
    shell:
        "{rules.clean_pre.shellcmd} && {rules.clean_conversions.shellcmd}"
