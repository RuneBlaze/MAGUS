# a version of mafft that is mafft-linsi and mafft add-frag in one
import argparse
from cProfile import run
from Bio import AlignIO, SeqIO
from tools import external_tools
from os.path import dirname, abspath
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", help="Input fasta file", required=True)
parser.add_argument("--output", "-o", help="Output fasta file", required=True)
parser.add_argument("--threads", "-t", type=int, help="# of threads", required=True)
parser.add_argument("-w", help="working directory")
parser.add_argument("--length", "-l", type=int, help="target length of sequences", required=True)
args = parser.parse_args()

if not args.w:
    args.w = dirname(args.output)

WORKING_DIR = args.w
NUM_THREADS = args.threads

def runMafftLinsi(input_p, output_p):
    external_tools.runMafft(input_p, None, WORKING_DIR, output_p, NUM_THREADS).run()

def runMafftFragAdd(ref_p, frag_p, output_p):
    external_tools.runMafftAdd(ref_p, frag_p, WORKING_DIR, output_p, NUM_THREADS, True).run()

def align_all(input_p, output_p, l):
    with open(input_p, 'r') as fh:
        alignment = list(SeqIO.parse(fh, 'fasta'))
    fulllength = []
    fragments = []
    for record in alignment:
        seq_l = len(record.seq.ungap("-"))
        if abs(seq_l - l) < 0.25 * l:
            fulllength.append(record)
        else:
            fragments.append(record)
    print("# of full length sequences for fafft: ", len(fulllength))
    fullPath = output_p + ".fl.fa"
    fragPath = output_p + ".nf.fa"
    SeqIO.write(fulllength, fullPath, "fasta")
    SeqIO.write(fragments, fragPath, "fasta")
    if fulllength and fragments:
        fullPathAligned = output_p + ".fl.aln"
        runMafftLinsi(abspath(fullPath), abspath(fullPathAligned))
        runMafftFragAdd(abspath(fullPathAligned), abspath(fragPath), abspath(output_p))
    elif fulllength:
        runMafftLinsi(abspath(fullPath), abspath(output_p))
    elif fragments:
        runMafftLinsi(abspath(fragPath), abspath(output_p))
    else:
        assert None

if __name__ == '__main__':
    align_all(args.input, args.output, args.length)