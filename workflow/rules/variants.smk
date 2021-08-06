
rule extract_links:
    output:
        pq=temp(paths.variants_batch.links_table),
    input:
        bam=paths.mapping.coord_sorted_bam,
        refgenome=paths.refgenome.fasta_unzipped,
        pore_c_table=paths.align_table.pore_c,
    params:
        vcf=lookup_value("vcf_path", mapping_df),
        dask_settings=config["software"]["dask"]["settings"],
    conda:
        PORE_C_CONDA_FILE
    threads: 1
    shell:
        "pore_c {params.dask_settings} --dask-num-workers {threads} "
        "variants extract-variant-info "
        "{input.bam} {params.vcf} {input.refgenome} {output.pq} "
        " --pore-c-table {input.pore_c_table} "


rule create_link_fofn:
    output:
        paths.variants.links_fofn,
    input:
        expand_basecall_batches(paths.variants_batch.links_table),
    run:
        with open(output[0], "w") as fh:
            fh.write("{}\n".format("\n".join(input)))


rule merge_links:
    input:
        fofn=paths.variants.links_fofn,
    output:
        pq=paths.variants.links_table,
    params:
        dask_settings=config["software"]["dask"]["settings"],
    conda:
        PORE_C_CONDA_FILE
    threads: 10
    run:
        raise NotImplementedError
        # "pore_c {params.dask_settings} --dask-num-workers {threads} "

        # "variants aggregate-links  "
        # "{input.fofn} {output.pq} --fofn "

