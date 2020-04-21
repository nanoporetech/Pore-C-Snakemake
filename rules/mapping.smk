rule align_bwa:
    output:
        bam = paths.mapping.coord_sorted_bam,
        bai = paths.mapping.coord_sorted_bai
    input:
        fastq = paths.basecall.fastq,
        refgenome = paths.refgenome.fasta,
        bwa_index = paths.refgenome.bwt,
    params:
        cli_opts = config['software']['bwa']['cli_opts'],
        memory = config['software']['sort']['memory_per_thread'],
        sort_threads = config['software']['sort']['threads']
    threads: config['software']['bwa']['threads']
    conda: "../envs/bwa.yml"
    log: to_log(paths.mapping.coord_sorted_bam) 
    benchmark: to_benchmark(paths.mapping.coord_sorted_bam) 
    shell:
        "( bwa {params.cli_opts} -t {threads} "
        "{input.refgenome}  {input.fastq} "
        " | python {workflow.basedir}/scripts/reformat_bam.py "   # TODO: replace this with pore_c command
        " | samtools sort -O bam -m {params.memory} -@ {params.sort_threads} -o {output.bam} -) 2>{log} ;"
        " samtools index {output.bam} 2>{log} "


rule haplotag:
    output:
        bam = paths.mapping.haplotagged_bam
    input:
        bam = paths.mapping.coord_sorted_bam,
        bai = paths.mapping.coord_sorted_bai,
        vcf = lookup_value('vcf_path', basecall_df),
        refgenome = paths.refgenome.fasta_unzipped
    conda: "../envs/whatshap.yml"
    log: to_log(paths.mapping.haplotagged_bam)
    benchmark: to_benchmark(paths.mapping.haplotagged_bam)
    shell:
        "whatshap haplotag -o {output} {input.vcf} {input.bam} --ignore-read-groups -r {input.refgenome} 2>{log}"


def is_phased(wildcards):
    vcf = lookup_value('vcf_path', basecall_df)(wildcards).strip()
    if vcf:
        return True
    else:
        return False

def choose_input_bam(wildcards):
    if is_phased(wildcards):
        return  paths.mapping.haplotagged_bam
    else:
        return  paths.mapping.coord_sorted_bam

def porec_phased_flag(wildcards):
    if is_phased(wildcards):
        return '--phased '
    else:
        return ' '


rule create_alignment_table:
    input:
        bam = choose_input_bam,
    output:
        pq = directory(paths.align_table.alignment)
    params:
        chunksize = 10000
    log: to_log(paths.align_table.alignment)
    benchmark: to_benchmark(paths.align_table.alignment)
    threads: config['software']['pore_c']['create_alignment_table']['threads']
    conda: PORE_C_CONDA_FILE
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        "alignments create-table {input.bam} {output.pq} --chunksize {params.chunksize} 2>{log}"

rule assign_fragments:
    input:
        align_table = paths.align_table.alignment,
        fragments_table = paths.virtual_digest.fragments
    output:
        pq = directory(paths.align_table.pore_c)
    log: to_log(paths.align_table.pore_c)
    benchmark: to_benchmark(paths.align_table.pore_c)
    threads: config['software']['pore_c']['create_alignment_table']['threads']
    conda: PORE_C_CONDA_FILE
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        "alignments assign-fragments {input.align_table} {input.fragments_table} {output} 2>{log}"


rule to_contacts:
    input:
        align_table = paths.align_table.pore_c,
    output:
        contacts = directory(paths.contacts.contacts),
        concatemers = directory(paths.contacts.concatemers)
    params:
        prefix = to_prefix(paths.contacts.contacts)
    log: to_log(paths.contacts.contacts)
    benchmark: to_benchmark(paths.contacts.contacts)
    conda: PORE_C_CONDA_FILE
    threads: 10
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        "alignments to-contacts {input} {output.contacts} {output.concatemers} 2>{log}"
