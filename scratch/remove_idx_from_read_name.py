from pathlib import Path
from pysam import AlignmentFile
from argparse import ArgumentParser

def remove_idx_from_read_names(
    input_bam: Path
):
    """ Replace READNAME:ALIGN_IDX with just READNAME

    Originally created because WhatsHap requires unique read names.
    """

    infile = AlignmentFile(input_bam, "rb")
    stdout = AlignmentFile("-", "wb", template=infile)
    align_iter = infile.fetch(until_eof=True)

    for read in align_iter:
        readname = read.query_name.split(":")[0]
        read.query_name = readname
        stdout.write(read)

    stdout.close()
    infile.close()

parser = ArgumentParser()
parser.add_argument("-i", "--input-bam", dest="input", help="bam to remove idx from")

args = parser.parse_args()

remove_idx_from_read_names(input_bam=args.input)

