rule create_base_matrix:
    input:
        pairs=paths.pairs.catalog,
    output:
        paths.matrix.catalog,
        paths.matrix.coo
    params:
        prefix=to_prefix(paths.matrix.catalog),
        resolution=config['matrix_resolutions']['base']
    benchmark: to_benchmark(paths.matrix.catalog)
    log: to_log(paths.matrix.catalog)
    threads: 1
    conda: PORE_C_CONDA_FILE
    shell:
        "pore_c pairs to-matrix {input} {params.prefix} -r {params.resolution} -n {threads} 2>{log}"


rule create_base_cooler_file:
    input:
        coo=paths.matrix.coo,
        chromsizes=paths.refgenome.chromsizes,
    output:
        paths.matrix.cool
    params:
        resolution=config['matrix_resolutions']['base'],
        refgenome_id=config['refgenome']['refgenome_id']
    benchmark: to_benchmark(paths.matrix.cool)
    log: to_log(paths.matrix.cool)
    conda: "../envs/cooler.yml"
    shell:
        "cooler load --assembly {params.refgenome_id} --input-copy-status unique -f coo {input.chromsizes}:{params.resolution}  {input.coo} {output} 2> {log}"


rule create_mcool_file:
    input:
        paths.matrix.cool
    output:
        paths.matrix.mcool_counts
    params:
        resolutions=",".join(map(str, config['matrix_resolutions']['zoomify']))
    benchmark: to_benchmark(paths.matrix.mcool_counts)
    log: to_log(paths.matrix.mcool_counts)
    conda: "../envs/cooler.yml"
    threads: 1
    shell:
        "cooler zoomify -n {threads} -r {params.resolutions} -o {output} {input} 2>{log}" 


