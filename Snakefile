import pandas as pd
import yaml
from pathlib import Path
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.5.4")

### Validation of schemas ###
#validate('config.yaml', schema="schemas/config.schema.yaml")
##### load config and sample sheets ##
configfile: "config.yaml"

wildcard_constraints: 
    enzyme="[^_]+",
    run_id="[^_]+",
    batch_id="batch\d+"

PORE_C_CONDA_FILE = "../envs/pore_c_dev.yml"

#validate(config['basecalls'], schema="schemas/basecall.schema.yaml")

basecall_df = (
        pd.read_table(config['basecalls'], comment='#')
        .set_index(['run_id', 'enzyme'], drop=False)
)

##### load rules #####
include: "rules/common.smk"  # python helper functions

outdir = Path(config['output_dir'])
paths = create_path_accessor(outdir)



include: "rules/refgenome.smk"  # prepare the reference genome
include: "rules/reads.smk"  # import fastqs
include: "rules/mapping.smk"  # map and process resulting alignments
include: "rules/exports.smk" # export to alternative formats

##### output paths #####



rule all:
    input:
        refgenome=paths.refgenome.catalog,
        contacts=expand_rows(paths.merged_contacts.concatemers, basecall_df)



rule cooler:
    input: expand_rows(paths.matrix.mcool, basecall_df)
    
rule haplotyped_cools:
    input: expand_rows(paths.matrix.haplotyped_cools, basecall_df)

rule pairs:
    input: expand_rows(paths.pairs.index, basecall_df)

rule salsa:
    input: expand_rows(paths.assembly.salsa_bed, basecall_df)

rule juicebox:
    input: expand_rows(paths.juicebox.hic, basecall_df)


rule test:
    input: 
        expand_rows(paths.matrix.mcool, basecall_df),
        expand_rows(paths.matrix.haplotyped_cools, basecall_df),
        expand_rows(paths.pairs.index, basecall_df),
        expand_rows(paths.assembly.salsa_bed, basecall_df),
        expand_rows(paths.juicebox.hic, basecall_df)
