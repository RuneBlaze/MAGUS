from Bio import AlignIO, SeqIO
import argparse
import os
from tools.external_tools import runFastTree, runOldRmEvaluate, runRaxmlEvaluate, runEpaNg, runGappaGraft, runPplacer


def estimate_tree(input, output):
    with open(input, 'r') as fh:
        alignment = AlignIO.read(fh, 'fasta')
    seqLengths = [len(record.seq.ungap("-")) for record in alignment]
    seqLengths.sort()
    topQuartile = seqLengths[int(0.75*(len(seqLengths)-1))]

    fulllength = []
    notfulllength = []

    for r in alignment:
        l = len(r.seq.ungap("-"))
        if abs(l - topQuartile) < 0.25 * topQuartile:
            fulllength.append(r)
        else:
            notfulllength.append(r)

    alignDir = os.path.abspath(output + "_env")
    def consumeTask(t):
        t.taskArgs["workingDir"] = alignDir
        t.run()

    fullPath = output + ".full"
    notFullPath = output + ".nf"
    if not os.path.exists(alignDir):
        os.mkdir(alignDir)
    SeqIO.write(fulllength, fullPath, "fasta")
    SeqIO.write(notfulllength, notFullPath, "fasta")
    fp = os.path.abspath(fullPath)
    consumeTask(runFastTree(os.path.abspath(fullPath), alignDir, os.path.abspath(fullPath + ".tre")))
    consumeTask(runOldRmEvaluate(os.path.abspath(fullPath), alignDir, fp + ".tre", os.path.abspath(fullPath + ".tre.eval")))
    consumeTask(runPplacer(os.path.abspath(fullPath), fp + ".tre.eval", os.path.abspath(notFullPath), alignDir, fp + ".jplace"))
    consumeTask(runGappaGraft(fp + ".jplace", alignDir, os.path.abspath(output)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
    description='A baseline method for estimating a good tree for a given alignment assuming sequence length heterogeneity')
    parser.add_argument('-i', '--input', help='Input alignment file')
    parser.add_argument('-o', '--output', help='Output tree file')

    args = parser.parse_args()
    estimate_tree(args.input, args.output)