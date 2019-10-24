rule coordinate_sort:
    output:
        bam = paths.phase.coordinate_sorted_bam,
        index = paths.phase.coordinate_sorted_bai
    input:
        bam = paths.mapping.read_sorted_bam
    params:
        memory = config['software']['sort']['memory_per_thread']
    threads: config['software']['sort']['threads']
    log: to_log(paths.phase.coordinate_sorted_bam)
    benchmark: to_benchmark(paths.phase.coordinate_sorted_bam)
    shell:
        "({ACTIVATE_POREC} python scratch/add_idx_to_read_name.py -i {input.bam} | "
        "samtools sort -O bam -m {params.memory} -@ {threads} -o {output.bam}; "
        "samtools index {output.bam}; ) 2>{log}"

rule unzip_reference_fasta:
    output: 
        fa = paths.refgenome.fasta_unzipped,
        fai = paths.refgenome.fasta_unzipped_fai
    input: gz = paths.refgenome.fasta 
    shell:
        "{ACTIVATE_POREC} gunzip -c {input.gz} > {output.fa};"
        "samtools faidx {output.fa};" 

rule haplotag:
    output:
        bam = paths.phase.haplotagged_bam
    input:
        bam = paths.phase.coordinate_sorted_bam,
        bai = paths.phase.coordinate_sorted_bai,
        vcf = config['sample']['phased_vcf_path'],
        refgenome = paths.refgenome.fasta_unzipped
    conda: "../envs/whatshap.yml"
    log: to_log(paths.phase.haplotagged_bam)
    benchmark: to_benchmark(paths.phase.haplotagged_bam)
    shell:
        "whatshap haplotag -o {output} {input.vcf} {input.bam} --ignore-read-groups -r {input.refgenome} 2>{log}"

rule re_read_sort:
    output:
        bam = paths.phase.read_sorted_haplotagged_bam
    input:
        bam = paths.phase.haplotagged_bam
    params:
        memory = config['software']['sort']['memory_per_thread']
    threads: config['software']['sort']['threads']
    log: to_log(paths.phase.read_sorted_haplotagged_bam)
    benchmark: to_benchmark(paths.phase.read_sorted_haplotagged_bam)
    shell:
        "({ACTIVATE_POREC} python scratch/remove_idx_from_read_name.py -i {input.bam} | "
        "samtools sort -O bam -n -m {params.memory} -@ {threads} -o {output.bam} )2>{log} "
