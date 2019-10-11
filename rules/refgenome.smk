rule add_refgenome:
    input: config['refgenome']['refgenome_path']
    output:
        catalog=paths.refgenome.catalog,
        chrom_metadata=paths.refgenome.chrom_metadata,
        chrom_sizes=paths.refgenome.chromsizes,
        fasta=paths.refgenome.fasta,
        fai=paths.refgenome.fai
    params:
        prefix=to_prefix(paths.refgenome.catalog),
        genome_id=config['refgenome']['refgenome_id']
    log: to_log(paths.refgenome.catalog)
    benchmark: to_benchmark(paths.refgenome.catalog)
    shell:
        "{ACTIVATE_POREC} pore_c refgenome catalog {input} {params.prefix} --genome-id {params.genome_id} 2> {log}"

rule virtual_digest:
    input:
        paths.refgenome.catalog
    output:
        paths.virtual_digest.catalog,
        paths.virtual_digest.fragments,
        paths.virtual_digest.digest_stats
    params:
        prefix=to_prefix(paths.virtual_digest.catalog)
    benchmark: to_benchmark(paths.virtual_digest.catalog)
    log: to_log(paths.virtual_digest.catalog)
    threads: 10
    shell:
        "{ACTIVATE_POREC} pore_c refgenome virtual-digest {input} enzyme:{wildcards.enzyme} {params.prefix} -n {threads} 2> {log}"

rule bwa_index_refgenome:
    input:
        paths.refgenome.fasta
    output:
        paths.refgenome.bwt
    conda: "../envs/bwa.yml"
    log: to_log(paths.refgenome.bwt)
    benchmark: to_benchmark(paths.refgenome.bwt)
    shell:
        "bwa index {input} 2>{log}"


