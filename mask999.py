# compress alignment. Removing columns that are insrtion only or 99.9% gaps.
import argparse
from Bio import AlignIO
from Bio.Seq import Seq
import numpy as np

def compress_alignment(input, output):
    A = AlignIO.read(open(input), "fasta")
    removable_idx = [] # removable column idx
    N = A.get_alignment_length()
    for i in range(N):
        row = A[:, i]
        if row.count('-') == len(row) - 1:
            removable_idx.append(i)
        elif row.count('-') >= len(row) * 0.999:
            removable_idx.append(i)
    for record in A:
        stringview = np.array(list(record.seq))
        stringview = np.delete(stringview, removable_idx)
        record.seq = Seq(''.join(list(stringview)))
    with open(output, "w+") as fh:
        AlignIO.write(A, fh, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Remove columns that are insertion only or 99.9% gaps, a la PASTA.')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)
    # assuming FASTA files
    args = parser.parse_args()
    compress_alignment(args.input, args.output)