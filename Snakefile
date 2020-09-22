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


BASE_DIR = Path(workflow.basedir)

### Validation of schemas ###
##### load config and sample sheets ##
configfile: BASE_DIR / "config/config.yaml"

##### load rules #####
include: "rules/common.smk" # python helper functions


basecall_df, reference_df, mapping_df = create_config_dataframes()
paths = create_path_accessor()


include: "rules/refgenome.smk" # prepare the reference genome


include: "rules/reads.smk" # import fastqs


include: "rules/mapping.smk" # map and process resulting alignments


include: "rules/exports.smk" # export to alternative formats


##### output paths #####


rule all:
    input:
        basecalls=expand_rows(paths.basecall.catalog, basecall_df),
        refgenome=expand_rows(paths.refgenome.bwt, reference_df),
        contacts=expand_rows(paths.merged_contacts.concatemers, mapping_df)


rule cooler:
    input:
        expand_rows(paths.matrix.mcool, mapping_df),


rule pairs:
    input:
        expand_rows(paths.pairs.index, mapping_df),


rule salsa:
    input:
        expand_rows(paths.assembly.salsa_bed, mapping_df),


rule juicer:
    input:
        expand_rows(paths.juicebox.hic, mapping_df),

rule mnd:
    input:
        expand_rows(paths.juicebox.mnd, mapping_df),



rule test:
    input:
        expand_rows(paths.merged_contacts.concatemers, mapping_df),
        expand_rows(paths.matrix.mcool, mapping_df),
        #expand_rows(paths.matrix.haplotyped_cools, mapping_df),
        expand_rows(paths.pairs.index, mapping_df),
        expand_rows(paths.assembly.salsa_bed, mapping_df),
        expand_rows(paths.juicebox.hic, mapping_df),
        expand_rows(paths.juicebox.mnd, mapping_df),
