import argparse
import os
import sys
import inspect
import subprocess
import glob
from glob import glob
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from tasks.remote_tasks import remote_mafft_linsi
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, help='Path to input unaligned sequences', required=True)
parser.add_argument('-s', '--support', type=int, help='Support size', required=True)
parser.add_argument('-o', '--output', type=str, help='Path to output alignment', required=True)
args = parser.parse_args()

INPUT_PATH = args.input
OUTPUT_PATH = args.output
ENV_PATH = os.path.relpath(args.input + '_env')
SUBTREE_PATH = INPUT_PATH + '.tre'
SUPPORT_SIZE = max(args.support, 100)
assert os.path.isfile(SUBTREE_PATH)
subprocess.run(['gcm137', 'slice', '-g', f'10x{SUPPORT_SIZE}', '-s', '100', '-i', INPUT_PATH, '-t', SUBTREE_PATH, '-o', ENV_PATH])
receipts = []
for unaligned_path in glob(os.path.join(ENV_PATH, "*", "*.unaln.fa")):
    print(f"Aligning: {unaligned_path}")
    d, bn = os.path.split(unaligned_path)
    aligned_path = os.path.join(d, bn.replace('.unaln.fa', '.aln.fa'))
    r = remote_mafft_linsi(os.path.abspath(unaligned_path), os.path.abspath(aligned_path))
    receipts.append(r)
for r in receipts:
    t = r.get(blocking=True)
    print(f"Alignment complete: {t.outputFile}")
aligned_constraints = glob(os.path.join(ENV_PATH, "constraints", "*.aln.fa"))
aligned_glues = glob(os.path.join(ENV_PATH, "glues", "*.aln.fa"))
subprocess.run(['gcm137', 'merge', '-i', *aligned_constraints, '-g', *aligned_glues, '-o', OUTPUT_PATH])