from snakemake.shell import shell
import sys

snakemake = snakemake  # noqa

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "a")

def write_empty_results():
    with open(snakemake.output[0], "w") as fh:
        fh.write("#readname\thaplotype\tphaseset\tchromosome\n")

# use a parameter to skip phasing if needed (makes snakemake pipeline easier to have empty file)
if not snakemake.params.is_phased:
    sys.stderr.write("Short-circuiting phasing information\n")
    write_empty_results()

    sys.exit(0)


# check if bam has no alignments
alignment_count = 0
for l in shell("samtools idxstats {snakemake.input.bam}", iterable=True):
    fields = l.split("\t")
    if fields[0] != "*":
        alignment_count += int(fields[-2])

if alignment_count == 0:
    sys.stderr.write(f"No mapped reads in bam {snakemake.input.bam}, writing empty phase results")
    write_empty_results()
    sys.exit(0)

shell(
    "whatshap haplotag  {snakemake.input.vcf} {snakemake.input.bam} "
    "--ignore-read-groups -r {snakemake.input.refgenome} "
    "--output-haplotag-list {snakemake.output}  "
    " -o /dev/null "
    " {log}"
)

