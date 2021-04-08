rule import_f5c_binary:
    output:
        gpu_binary=config["software"]["f5c"]["gpu_binary"],
        cpu_binary=config["software"]["f5c"]["cpu_binary"],
    params:
        f5curl=config["software"]["f5c"]["tools_url"],
        version=config["software"]["f5c"]["version"],
    shell:
        """
        curl -L {params.f5curl} \
            | tar -zxf - ;
        mv f5c-{params.version} bin/f5c ;
        """


rule filter_bam:
    # Returns a filtered bam, containing only passing alignments.
    input:
        bam=paths.mapping.coord_sorted_bam,
        porec_table=paths.align_table.pore_c,
    output:
        filtered_bam=paths.mapping.filtered_bam,
        filtered_bai=paths.mapping.filtered_bai,
    params:
        script="scripts/filter_bam.py",
    benchmark:
        to_benchmark(paths.mapping.filtered_bam)
    log:
        to_log(paths.mapping.filtered_bam),
    conda:
        "../envs/python_script_utils.yml"
    shell:
        """
        ( python {params.script} -i {input.bam} --porec-table {input.porec_table} -o {output.filtered_bam} ;
        samtools index {output.filtered_bam} ; ) 2>{log}
        """


rule f5c_index:
    # TODO: f5c_index creates a .readdb containing a mapping for all reads,
    # not just reads that are in the fastq. Since we are working with batched
    # fqs but not f5s this adds some extra time. Indexing 1 batch takes ~1 minute.
    # TODO: Hasindu recently released a faster way to index when lacking a
    # sequencing summary file, which isn't implemented here.
    input:
        f5c=config["software"]["f5c"]["cpu_binary"],
        fast5=paths.fast5.fast5,
        fastq=paths.basecall.fastq,
        summary=paths.fast5.seq_summary_pf,
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
        {input.f5c} index -d {input.fast5} -s {input.summary} {input.fastq} 2>{log}
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
        f5c=config["software"]["f5c"]["binary"],
        bam=paths.mapping.filtered_bam,
        index=paths.basecall.f5c_index,
        fastq=paths.basecall.fastq,
        reference=paths.refgenome.fasta,
    output:
        per_read=paths.methylation.per_read_llr,
    resources:
        gpu=config["software"]["f5c"]["gpus"],
    params:
        settings=config["software"]["f5c"]["cli_opts"],
    benchmark:
        to_benchmark(paths.methylation.per_read_llr)
    log:
        to_log(paths.methylation.per_read_llr),
    threads: config["software"]["f5c"]["threads"]
    shell:
        """
        {input.f5c} call-methylation -t {threads} \
            {params.settings} -B 6.0M -K 550 --iop 20 \
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
        f5c=config["software"]["f5c"]["cpu_binary"],
        per_read=paths.methylation.combined_read_llr,
    output:
        split_pos=paths.methylation.per_locus_methylation,
    benchmark:
        to_benchmark(paths.methylation.per_locus_methylation)
    log:
        to_log(paths.methylation.per_locus_methylation),
    shell:
        """
        {input.f5c} meth-freq -s -i {input.per_read} > {output.split_pos} 2>{log}
        """
