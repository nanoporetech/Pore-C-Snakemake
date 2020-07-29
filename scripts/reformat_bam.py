import sys
import pysam


def reformat_bam(infh, outfh):

    infile = pysam.AlignmentFile(infh, "r")
    outfile = pysam.AlignmentFile(outfh, "w", template=infile)
    read_indices = {}
    for align_idx, align in enumerate(infile.fetch(until_eof=True)):
        read_id = align.query_name
        read_idx = read_indices.get(read_id, None)
        if read_idx is None:
            read_idx = len(read_indices)
            read_indices[read_id] = read_idx
        # align.set_tag(tag='BX', value=align.query_name, value_type="Z")
        align.query_name = f"{read_id}:{read_idx}:{align_idx}"
        outfile.write(align)


if __name__ == "__main__":
    reformat_bam(sys.stdin, sys.stdout)
