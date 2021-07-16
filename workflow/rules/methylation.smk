rule import_f5c_binary:
    output:
        gpu_binary=paths.binaries.f5c_gpu,
        cpu_binary=paths.binaries.f5c_cpu,
    params:
        f5curl=config["software"]["f5c"]["tools_url"],
        version=config["software"]["f5c"]["version"],
        bindir=paths.binaries.f5c_cpu.rsplit("/", 1)[0],
    shell:
        """
        mkdir -p {params.bindir} ;
        curl -L {params.f5curl} \
            | tar -zxf - ;
        mv f5c-{params.version}/f5c_x86_64_linux* {params.bindir};
        """


rule filter_bam:
    # Returns a filtered bam, containing only passing alignments.
    input:
        bam=paths.mapping.coord_sorted_bam,
        pore_c_table=paths.align_table.pore_c,
    output:
        filtered_bam=paths.mapping.filtered_bam,
        filtered_bai=paths.mapping.filtered_bai,
    params:
        dask_settings=config["software"]["dask"]["settings"],
    benchmark:
        to_benchmark(paths.mapping.filtered_bam)
    log:
        to_log(paths.mapping.filtered_bam),
    conda:
        PORE_C_CONDA_FILE
    shell:
        "( pore_c {params.dask_settings} --dask-num-workers {threads} "
        "alignments filter-bam {input.bam} {input.pore_c_table} {output.filtered_bam} "
        " --clean-read-name ;"
        " samtools index {output.filtered_bam} ; ) 2>{log} "


f5c_mode = config["software"]["f5c"]["run_mode"]
f5c_config = config["software"]["f5c"]["settings"][f5c_mode]


rule f5c_index:
    # TODO: f5c_index creates a .readdb containing a mapping for all reads,
    # not just reads that are in the fastq. Since we are working with batched
    # fqs but not f5s this adds some extra time. Indexing 1 batch takes ~1 minute.
    # TODO: Hasindu recently released a faster way to index when lacking a
    # sequencing summary file, which isn't implemented here.
    input:
        fast5=paths.fast5.fast5,
        fastq=paths.basecall.fastq,
        summary=paths.fast5.seq_summary_pf,
        binary=paths.binaries[f"f5c_{f5c_mode}"],
    output:
        index=paths.basecall.f5c_index,
        fai=paths.basecall.f5c_fai,
        gzi=paths.basecall.f5c_gzi,
        readdb=paths.basecall.f5c_readdb,
    benchmark:
        to_benchmark(paths.basecall.f5c_index)
    log:
        to_log(paths.basecall.f5c_index),
    shell:
        """
        {input.binary} index -d {input.fast5} -s {input.summary} {input.fastq} 2>{log}
        """


rule f5c_index_all_files:
    input:
        expand_basecall_batches(paths.basecall.f5c_index),
    output:
        paths.basecall.f5c_index_all,
    shell:
        "touch {output}"


rule f5c_call_methylation:
    # CPU or GPU is inelegantly determined by commenting out lines in config.yaml
    input:
        bam=paths.mapping.filtered_bam,
        index=paths.basecall.f5c_index,
        fastq=paths.basecall.fastq,
        reference=paths.refgenome.fasta,
        binary=paths.binaries[f"f5c_{f5c_mode}"],
    output:
        per_read=paths.methylation.per_read_llr,
    resources:
        gpu=f5c_config["gpus"],
    params:
        cli_opts=f5c_config["cli_opts"],
    benchmark:
        to_benchmark(paths.methylation.per_read_llr)
    log:
        to_log(paths.methylation.per_read_llr),
    threads: f5c_config["threads"]
    shell:
        """
        {input.binary} call-methylation -t {threads} \
            {params.cli_opts} -B 6.0M -K 550 --iop 20 \
            -b {input.bam} -g {input.reference} -r {input.fastq} \
            > {output.per_read} 2>{log}
        """


rule f5c_aggregate:
    input:
        expand_basecall_batches(paths.methylation.per_read_llr),
    output:
        paths.methylation.combined_read_llr,
    run:
        shell("head -n 1 {input[0]} > {output} ")
        for f in input:
            shell("tail -n+2 {f}  >> {output}")


rule f5c_meth_freq:
    # Takes per-read, per-position calls and summarizes them to get per-position methylation frequencies.
    input:
        per_read=paths.methylation.combined_read_llr,
        binary=paths.binaries[f"f5c_{f5c_mode}"],
    output:
        split_pos=paths.methylation.per_locus_methylation,
    params:
        f5c=f5c_config["binary"],
    benchmark:
        to_benchmark(paths.methylation.per_locus_methylation)
    log:
        to_log(paths.methylation.per_locus_methylation),
    shell:
        """
        {input.binary} meth-freq -s -i {input.per_read} > {output.split_pos} 2>{log}
        """
