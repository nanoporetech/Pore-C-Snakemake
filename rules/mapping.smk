rule align_bwa:
    output:
        bam = paths.mapping.read_sorted_bam
    input:
        fastq = paths.basecall.fastq,
        refgenome = paths.refgenome.fasta,
        bwa_index = paths.refgenome.bwt,
    params:
        cli_opts = config['software']['bwa']['cli_opts']
    threads: config['software']['bwa']['threads']
    conda: "../envs/bwa.yml"
    log: to_log(paths.mapping.read_sorted_bam) 
    benchmark: to_benchmark(paths.mapping.read_sorted_bam) 
    shell:
        "( bwa {params.cli_opts} -t {threads} "
        "{input.refgenome}  {input.fastq} "
        " | samtools sort -O bam -m 4G -n -o {output} -) 2>{log}"


rule create_alignment_table:
    input:
        bam = paths.mapping.read_sorted_bam,
        virtual_digest = paths.virtual_digest.catalog
    output:
        paths.align_table.catalog,
        paths.align_table.alignment,
        paths.align_table.read,
        paths.align_table.overlap,
        paths.align_table.alignment_summary,
        paths.align_table.read_summary
    params:
        prefix=to_prefix(paths.align_table.catalog)
    log: to_log(paths.align_table.catalog)
    benchmark: to_benchmark(paths.align_table.catalog)
    threads: config['software']['pore_c']['create_alignment_table']['threads']
    shell:
        "{ACTIVATE_POREC} pore_c alignments parse {input.bam} {input.virtual_digest} {params.prefix} -n {threads} 2>{log}"
