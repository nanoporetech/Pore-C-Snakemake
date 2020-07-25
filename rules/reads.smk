checkpoint import_basecalls:
    output:
        paths.basecall.catalog,
        paths.basecall.summary,
        paths.basecall.read_metadata
    params:
        fname = lookup_value('fastq_path', basecall_df),
        prefix = to_prefix(paths.basecall.catalog),
    log: to_log(paths.basecall.catalog)
    benchmark: to_benchmark(paths.basecall.catalog)
    conda: PORE_C_CONDA_FILE
    threads: 1
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        "reads prepare {params.fname} {params.prefix} --batch-size {config[reads_per_batch]} 2> {log}"



