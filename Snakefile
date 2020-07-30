import pandas as pd
import yaml
from pathlib import Path
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.5.4")


wildcard_constraints:
    enzyme="[^_]+",
    run_id="[^_]+",
    batch_id="batch\d+",


### Validation of schemas ###
##### load config and sample sheets ##
configfile: "config.yaml"


if config["pore_c_version"] == "rel":
    PORE_C_CONDA_FILE = "../envs/pore_c_rel.yml"
else:
    assert config["pore_c_version"] == "dev"
    PORE_C_CONDA_FILE = "../envs/pore_c_dev.yml"


basecall_df = pd.read_table(config["basecalls"], comment="#").set_index(["run_id", "enzyme"], drop=False)


##### load rules #####
include: "rules/common.smk" # python helper functions


outdir = Path(config["output_dir"])
paths = create_path_accessor(outdir)


include: "rules/refgenome.smk" # prepare the reference genome


include: "rules/reads.smk" # import fastqs


include: "rules/mapping.smk" # map and process resulting alignments


include: "rules/exports.smk" # export to alternative formats


##### output paths #####


rule all:
    input:
        refgenome=paths.refgenome.catalog,
        contacts=expand_rows(paths.merged_contacts.concatemers, basecall_df),


rule cooler:
    input:
        expand_rows(paths.matrix.mcool, basecall_df),


rule haplotyped_cools:
    input:
        expand_rows(paths.matrix.haplotyped_cools, basecall_df),


rule pairs:
    input:
        expand_rows(paths.pairs.index, basecall_df),


rule salsa:
    input:
        expand_rows(paths.assembly.salsa_bed, basecall_df),


rule juicebox:
    input:
        expand_rows(paths.juicebox.hic, basecall_df),


rule test:
    input:
        expand_rows(paths.merged_contacts.concatemers, basecall_df),
        expand_rows(paths.matrix.mcool, basecall_df),
        expand_rows(paths.matrix.haplotyped_cools, basecall_df),
        expand_rows(paths.pairs.index, basecall_df),
        expand_rows(paths.assembly.salsa_bed, basecall_df),
        expand_rows(paths.juicebox.hic, basecall_df),
