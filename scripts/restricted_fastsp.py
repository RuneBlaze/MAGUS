import argparse
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
import numpy as np
import os

parser = argparse.ArgumentParser(description='Run FastSP on a restricted set of sequences.')
parser.add_argument('-r', '--ref', type=str, required=True, help="reference MSA file")
parser.add_argument('-e', '--est', type=str, required=True, help="estimated MSA file")

ref_ids = set()

args = parser.parse_args()
ref_alignment = AlignIO.read(args.ref, "fasta")
for record in ref_alignment:
    ref_ids.add(record.id)

# filter the estimated alignment
est_alignment = AlignIO.read(args.est, "fasta")
filtered = [record for record in est_alignment if record.id in ref_ids]
est_ids = set([record.id for record in est_alignment if record.id in ref_ids])
filtered_ref = [record for record in ref_alignment if record.id in est_ids]

def compress_alignment(input, output):
    A = AlignIO.read(open(input), "fasta")
    removable_idx = [] # removable column idx
    N = A.get_alignment_length()
    for i in range(N):
        row = A[:, i]
        if row.count('-') == len(row):
            removable_idx.append(i)
    for record in A:
        stringview = np.array(list(record.seq))
        stringview = np.delete(stringview, removable_idx)
        record.seq = Seq(''.join(list(stringview)))
    with open(output, "w+") as fh:
        AlignIO.write(A, fh, "fasta")

if len(filtered) <= 1:
  print("")
  exit()

# write out the estimated alignment
# restricted_ = args.est + ".restricted_"
restricted = args.est + ".restricted"
restricted_ref = args.est + ".rref"
SeqIO.write(filtered, restricted, "fasta")
SeqIO.write(filtered_ref, restricted_ref, "fasta")
compress_alignment(restricted_ref, restricted_ref + "_")
os.system("java -jar FastSP.jar -e %s -r %s" % (restricted, restricted_ref + "_"))
print(f"SeqNum {len(filtered)}")