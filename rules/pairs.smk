rule create_pairs_file:
    input:
        paths.align_table.catalog
    output:
        catalog=paths.pairs.catalog,
        pairs=paths.pairs.pairs,
    params:
        prefix=to_prefix(paths.pairs.catalog)
    benchmark: to_benchmark(paths.pairs.catalog)
    log: to_log(paths.pairs.catalog)
    threads: config['software']['pore_c']['create_pairs_file']['threads']
    shell:
        "{ACTIVATE_POREC} pore_c pairs from-alignment-table {input} {params.prefix} -n {threads} 2>{log}"
