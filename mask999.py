import argparse
from numba import jit
import numpy as np
from Bio import SeqIO

# fast alignment masker that uses numba
# and also takes into account 'X' characters for protein datasets
# (TODO: handle NT data)

@jit(nopython=True)
def count_gaps(seq, gapcounts):
    for i, c in enumerate(seq):
        if c == 45 or c == 88:
            gapcounts[i] += 1

def compress_alignment(input, output):
    tt = 0
    for i, record in enumerate(SeqIO.parse(input, 'fasta')):
        ba = bytearray(record.seq._data)
        if i == 0:
            C = len(ba)
            gapcounts = np.zeros(C, dtype=np.int32)
        count_gaps(ba, gapcounts)
        tt += 1
    threshold = tt * 0.999
    with open(output, 'w+') as fh:
        for i, record in enumerate(SeqIO.parse(input, 'fasta')):
            fh.write(">{}\n".format(record.id))
            for j, c in enumerate(record.seq._data):
                if gapcounts[j] < threshold:
                    fh.write(chr(c))
            fh.write("\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Remove columns that are insertion only or 99.9% gaps, a la PASTA.')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)
    # assuming FASTA files
    args = parser.parse_args()
    compress_alignment(args.input, args.output)