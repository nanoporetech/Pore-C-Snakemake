import pandas as pd
import yaml
from pathlib import Path
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
#min_version("5.3.0")


##### load config and sample sheets ##

configfile: "config.yaml"

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
include: "rules/pairs.smk"  # convert alignments to pairs files 
include: "rules/matrices.smk" # convert to .cool and .mcool files


##### output paths #####



rule all:
    input: 
        refgenome=paths.refgenome.catalog,
        matrix=expand_rows(paths.matrix.mcool_counts, basecall_df)
    
        #reference_fasta=API.refgenome.all('bwt'),
        #virtual_digest=API.virtual_digest.all('catalog'),
        #mapping=API.mapping.all("read_sorted_bam", query),
        #align_table=API.align_table.all("catalog", query),
        ##pairs=API.pairs.all("catalog", query),
        ##balanced_matrix=API.matrix.all("mcool_balanced", query),
        ##correlation=API.refdataset_compare_correlate.all('catalog', query),
        ##tads = API.tad_results.all('sentinel', query),
        ##merged_matrix_count = API.rungroup_matrix.all("mcool_counts"),
        #merged_matrix_balanced = API.rungroup_matrix.all("mcool_balanced"),
        ##pairs=API.pairs.all("catalog", query),
        #balanced_matrix=API.matrix.all("mcool_balanced", query),
        #rungroup_matrix_downsample = API.rungroup_matrix_downsample.all("cool")




