checkpoint import_basecalls:
    output:
        paths.basecall.catalog,
        paths.basecall.summary,
        paths.basecall.read_metadata,
    params:
        fname=lookup_value("fastq_path", basecall_df),
        prefix=to_prefix(paths.basecall.catalog),
    log:
        to_log(paths.basecall.catalog),
    benchmark:
        to_benchmark(paths.basecall.catalog)
    conda:
        PORE_C_CONDA_FILE
    threads: 1
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        "reads prepare {params.fname} {params.prefix} --max-read-length {config[max_read_length]} "
        " --batch-size {config[reads_per_batch]} 2> {log}"


rule import_fast5s:
    # Makes symlinks to the fast5 files, which can then be renamed for
    # for f5c index without modifying the original fast5s.
    output:
        directory(paths.fast5.fast5),
    params:
        fname=lookup_value("fast5_directory", basecall_df),
    run:
        expand_tilde = os.path.expanduser(params.fname)
        path_is_absolute = os.path.isabs(expand_tilde)
        # cp -rs requires an absolute path for the source
        if path_is_absolute:
            shell("cp -rs {params.fname} {output}")
        else:
            shell('cp -rs "$(pwd)/{params.fname}" {output}')
