
rule to_salsa_bed:
    output: paths.assembly.salsa_bed
    input: paths.align_table.catalog
    log: to_log(paths.assembly.salsa_bed)
    benchmark: to_benchmark(paths.assembly.salsa_bed)
    threads: config['software']['pore_c']['to_salsa_bed']['threads']
    shell:
        "{ACTIVATE_POREC} pore_c alignments to-salsa-bed {input} {output} -n {threads}"


rule create_hicRef:
    output: paths.juicebox.hicref
    input: paths.virtual_digest.catalog
    log: to_log(paths.juicebox.hicref)
    benchmark: to_log(paths.juicebox.hicref)
    shell:
        "{ACTIVATE_POREC} pore_c refgenome to-hicref {input} {output} 2>{log}"

rule create_hic_txt:
    output: paths.juicebox.hic_txt
    input:
        hicref=paths.juicebox.hicref, # dummy dependency, file needed for conversion to binary
        align_catalog=paths.align_table.catalog
    log: to_log(paths.juicebox.hic_txt)
    shell:
       "{ACTIVATE_POREC} pore_c alignments to-hic-txt {input.align_catalog} {output} 2>{log}"
