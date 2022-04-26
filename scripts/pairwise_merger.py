import argparse
import treeswift as ts
import os
parser = argparse.ArgumentParser(description='Pairwise merge sequences using GCM')
parser.add_argument('-i', '--input')
# parser.add_argument('-o', '--output')

args = parser.parse_args()

MST = ts.read_tree_newick(args.input)

basepath = os.path.dirname(os.path.dirname(args.input))

graphbasepath = os.path.join(basepath, "graph")
alnbasepath = os.path.join(basepath, "subalignments")

def nthsupp(i):
    return os.path.join(graphbasepath, f"backbone_{i}_mafft.txt")

def nthsubset(i):
    return os.path.join(alnbasepath, f"subalignment_subset_{i}.txt")

def nthoutput(u, v):
    a = min(u, v)
    b = max(u, v)
    return os.path.join(basepath, "merged", f"merged_{a}_{b}.txt")

def dfs_edges(tree):
    results = []
    start = next(tree.traverse_leaves())
    stack = [start]
    visited = set([start.label])
    reorder = lambda tup: (min(tup), max(tup))
    while stack:
        cur = stack.pop()
        if not cur.is_leaf():
            for child in cur.children:
                if child.label not in visited:
                    stack.append(child)
                    results.append(reorder((int(cur.label), int(child.label))))
                    visited.add(child.label)
        if not cur.is_root():
            if cur.parent.label not in visited:
                stack.append(cur.parent)
                if cur.label and cur.parent.label:
                    results.append(reorder((int(cur.label), int(cur.parent.label))))
                visited.add(cur.parent.label)
    return results

# print(dfs_edges(MST), len(dfs_edges(MST)))

nodes = set([])
edges = []
for n in MST.traverse_postorder():
    if n.is_root():
        continue
    if not n.label:
        assert False
    edges.append((int(n.label), int(n.parent.label)))
    nodes.add(int(n.label))
    nodes.add(int(n.parent.label))


import pathlib
pathlib.Path(os.path.join(basepath, "merged")).mkdir(parents=True, exist_ok=True) 
# os.mkdir(, exist_ok=True)

with open(os.path.join(basepath, "merged", "merger.txt"), "w+") as fh:
    for i, (u, v) in enumerate(edges):
        u_ss = nthsubset(u)
        v_ss = nthsubset(v)
        i_supp = nthsupp(i+1)
        fh.write(f"{u_ss} {v_ss} {i_supp} {nthoutput(u, v)}\n")


with open(os.path.join(basepath, "merged", "order.txt"), "w+") as fh:
    for u, v in dfs_edges(MST):
        fh.write(nthoutput(u, v) + "\n")

# for u, v in dfs_edges(MST):
#     print(nthoutput(u, v))