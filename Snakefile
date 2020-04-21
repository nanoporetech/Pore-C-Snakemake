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
#include: "rules/pairs.smk"  # convert alignments to pairs files
#include: "rules/matrices.smk" # convert to .cool and .mcool files
include: "rules/exports.smk" # export to alternative formats
#include: "rules/phase.smk" # optional: phase the alignments

##### output paths #####



rule all:
    input:
        refgenome=paths.refgenome.catalog,
        contacts=expand_rows(paths.contacts.contacts, basecall_df)
        #matrix=expand_rows(paths.matrix.mcool_counts, basecall_df)
#
#rule phased:
#    input:
#        haplotagged=expand_rows(paths.phase.read_sorted_haplotagged_bam, basecall_df)
#

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
