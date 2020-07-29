from pathlib import Path
from pysam import AlignmentFile
from argparse import ArgumentParser


def add_idx_to_read_name(input_bam: Path):
    """ Changes the readname to be READNAME:ALIGN_IDX to have 'unique' readnames

    WhatsHap requires unique read names.
    """

    infile = AlignmentFile(input_bam, "rb")
    stdout = AlignmentFile("-", "wb", template=infile)
    align_iter = infile.fetch(until_eof=True)

    i = 0
    for read in align_iter:
        read.query_name = read.query_name + ":" + str(i)
        stdout.write(read)
        i = i + 1

    stdout.close()
    infile.close()


parser = ArgumentParser()
parser.add_argument("-i", "--input-bam", dest="input", help="bam to add idx to")

args = parser.parse_args()

add_idx_to_read_name(input_bam=args.input)
