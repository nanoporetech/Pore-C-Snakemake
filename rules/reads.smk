rule import_basecalls:
    output:
        paths.basecall.catalog,
        paths.basecall.fastq,
        paths.basecall.summary,
        paths.basecall.read_metadata
    params:
        fname = lookup_value('fastq_path', basecall_df),
        prefix = to_prefix(paths.basecall.catalog)
    log: to_log(paths.basecall.catalog)
    benchmark: to_benchmark(paths.basecall.catalog)
    shell:
        "{ACTIVATE_POREC} pore_c reads catalog  {params.fname} {params.prefix} 2> {log}"
