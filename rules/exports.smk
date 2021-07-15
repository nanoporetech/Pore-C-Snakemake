rule to_cooler:
    output:
        paths.matrix.cool,
    input:
        contacts=paths.merged_contacts.contacts,
        chromsizes=paths.refgenome.chromsizes,
        fragments=paths.virtual_digest.fragments,
    params:
        prefix=to_prefix(paths.matrix.cool, 1),
    log:
        to_log(paths.matrix.cool),
    benchmark:
        to_benchmark(paths.matrix.cool)
    threads: config["software"]["pore_c"]["to_cooler"]["threads"]
    conda:
        PORE_C_CONDA_FILE
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        " contacts export {input.contacts} cooler {params.prefix} --fragment-table {input.fragments} --chromsizes {input.chromsizes} 2>{log} "


rule to_haplotyped_cooler:
    output:
        touch(paths.matrix.haplotyped_cools),
    input:
        contacts=paths.merged_contacts.contacts,
        chromsizes=paths.refgenome.chromsizes,
        fragments=paths.virtual_digest.fragments,
    params:
        prefix=to_prefix(paths.matrix.haplotyped_cools, 2),
    log:
        to_log(paths.matrix.haplotyped_cools),
    benchmark:
        to_benchmark(paths.matrix.haplotyped_cools)
    threads: config["software"]["pore_c"]["to_cooler"]["threads"]
    conda:
        PORE_C_CONDA_FILE
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        " contacts export {input.contacts} cooler {params.prefix} --by-haplotype --fragment-table {input.fragments} --chromsizes {input.chromsizes} 2>{log} "


rule create_mcool_file:
    input:
        paths.matrix.cool,
    output:
        paths.matrix.mcool,
    params:
        resolutions=",".join(map(str, config["matrix_resolutions"]["zoomify"])),
    benchmark:
        to_benchmark(paths.matrix.mcool)
    log:
        to_log(paths.matrix.mcool),
    conda:
        "../envs/cooler.yml"
    threads: 1
    shell:
        "cooler zoomify -n {threads} -r {params.resolutions} -o {output} {input} 2>{log}"


rule to_unsorted_pairs:
    output:
        paths.pairs.unsorted_pairs,
    input:
        contacts=paths.merged_contacts.contacts,
        chromsizes=paths.refgenome.chromsizes,
    params:
        prefix=to_prefix(paths.pairs.unsorted_pairs, 1),
    log:
        to_log(paths.pairs.unsorted_pairs),
    benchmark:
        to_benchmark(paths.pairs.unsorted_pairs)
    threads: config["software"]["pore_c"]["to_unsorted_pairs"]["threads"]
    conda:
        PORE_C_CONDA_FILE
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        " contacts export {input.contacts} pairs {params.prefix} --chromsizes {input.chromsizes} 2>{log} "


rule to_sorted_pairs:
    output:
        paths.pairs.pairs,
    input:
        paths.pairs.unsorted_pairs,
    params:
        uncomp_file=to_prefix(paths.pairs.pairs, 1),
    log:
        to_log(paths.pairs.pairs),
    benchmark:
        to_benchmark(paths.pairs.pairs)
    threads: config["software"]["pore_c"]["sort_pairs_file"]["threads"]
    conda:
        "../envs/pairtools.yml"
    shell:
        "(pairtools sort {input} --nproc {threads} > {params.uncomp_file} &&  bgzip {params.uncomp_file} -) 2>{log}"


rule index_pairs:
    output:
        paths.pairs.index,
    input:
        paths.pairs.pairs,
    log:
        to_log(paths.pairs.index),
    benchmark:
        to_benchmark(paths.pairs.index)
    conda:
        "../envs/pairix.yml"
    shell:
        "pairix {input} 2>{log}"


rule to_salsa_bed:
    output:
        paths.assembly.salsa_bed,
    input:
        contacts=paths.merged_contacts.contacts,
    params:
        prefix=to_prefix(paths.assembly.salsa_bed),
    log:
        to_log(paths.assembly.salsa_bed),
    benchmark:
        to_benchmark(paths.assembly.salsa_bed)
    threads: config["software"]["pore_c"]["to_salsa_bed"]["threads"]
    conda:
        PORE_C_CONDA_FILE
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        " contacts export {input.contacts} salsa_bed {params.prefix}  2>{log}"


rule download_juicer_tools:
    output:
        paths.juicebox.tools,
    params:
        url=config["software"]["juicer"]["tools_url"],
    conda:
        "../envs/juicer_tools.yml"
    log:
        to_log(paths.juicebox.tools),
    shell:
        "wget -O - {params.url} > {output} 2>{log}"


rule create_hicRef:
    output:
        paths.juicebox.hicref,
    input:
        paths.virtual_digest.fragments,
    log:
        to_log(paths.juicebox.hicref),
    benchmark:
        to_log(paths.juicebox.hicref)
    conda:
        PORE_C_CONDA_FILE
    threads: 1
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        "refgenome fragments-to-hicref {input} {output} 2>{log}"


rule create_hic:
    output:
        paths.juicebox.hic,
    input:
        hicref=paths.juicebox.hicref,  # dummy dependency, file needed for conversion to binary
        chromsizes=paths.refgenome.chromsizes,
        pairs=paths.pairs.pairs,
        tools=paths.juicebox.tools,
    log:
        to_log(paths.juicebox.hic),
    conda:
        "../envs/juicer_tools.yml"
    shell:
        "java -Xmx2g -jar {input.tools} pre {input.pairs} {output} {input.chromsizes} -f {input.hicref} &>{log}"


rule to_mnd:
    output:
        paths.juicebox.mnd,
    input:
        contacts=paths.merged_contacts.contacts,
        refgenome=paths.refgenome.fasta,
    params:
        prefix=to_prefix(paths.juicebox.mnd),
    log:
        to_log(paths.juicebox.mnd),
    benchmark:
        to_benchmark(paths.juicebox.mnd)
    threads: 1
    conda:
        PORE_C_CONDA_FILE
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        " contacts export {input.contacts} merged_no_dups {params.prefix} --reference-fasta {input.refgenome} 2>{log}"
