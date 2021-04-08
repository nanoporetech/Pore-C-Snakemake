import os
import glob
from pathlib import Path
from pysam import AlignmentFile
import pyarrow
from pyarrow import parquet as pq
import pandas as pd
from natsort import natsorted, ns
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--input-bam", dest="input",
    help="Unfiltered bam")
parser.add_argument("-o", "--output-bam", dest="output",
    help="Destination to write the filtered bam to")
parser.add_argument("-p", "--porec-table", dest="porec",
    help="Directory of fragment alignment tables containing best tilepaths")

args = parser.parse_args()

# def read_porec_parquet(path):
#     parquet = pd.read_table(path).to_pandas()
#     parquet.sort_values('align_idx', inplace=True)
#     parquet.reset_index(inplace=True)

# def read_all_parquets(directory):
#     all_paths = natsorted(glob.glob(os.path.join(directory, "*.parquet")))
#     generator = (pq.read_table(p) for p in all_paths)
#     parquets = pyarrow.concat_tables(generator).to_pandas()
#     parquets.sort_values('align_idx', inplace=True)
#     parquets.reset_index(inplace=True)
#     if all(parquets.index == parquets['align_idx']) == False:
#         raise(Exception("Index doesn't match align_idx ; are some parquets missing?"))
#     return(parquets)

def filter_bam(input_bam, output_bam, parquets):
    inbam = AlignmentFile(input_bam, "rb")
    outbam = AlignmentFile(output_bam, "wb", template=inbam)
    for align in inbam.fetch(until_eof=True):
        align_idx = int(align.query_name.rsplit(":")[2])
        filter_status = parquets.iloc[align_idx]['pass_filter']
        if filter_status == True:
            readname_only = align.query_name.split(":")[0]
            align.query_name = readname_only
            outbam.write(align)

if __name__ == '__main__':
    parquets = pq.read_table(args.porec).to_pandas()
    parquets.sort_values('align_idx', inplace=True)
    filter_bam(args.input, args.output, parquets)
