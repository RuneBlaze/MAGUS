from collections import defaultdict
def edges(t):
    seen = set()
    for n in t.traverse_postorder():
        if n.is_root():
            continue
        u = int(n.label)
        v = int(n.parent.label)
        r = min(u, v), max(u, v)
        if r not in seen:
            seen.add(r)
            yield r

def evolve(frontier, adjacencies):
    for u in frontier:
        for v in adjacencies[u]:
            if v not in frontier:
                yield frontier | frozenset([v])

def enumerate_subgraphs(t, size_limit):
    "t : Treeswift.Tree"

    adjacencies = defaultdict(set)
    # now stack contains the singleton nodes
    for u, v in edges(t):
        adjacencies[u].add(v)
        adjacencies[v].add(u)
    stack = []
    visited = set()
    for n in t.traverse_postorder():
        if n.label:
            singleton = frozenset([int(n.label)])
            stack.append(singleton)
            visited.add(singleton)
    while stack:
        cur = stack.pop()
        for next_frontier in evolve(cur, adjacencies):
            if next_frontier not in visited:
                visited.add(next_frontier)
                if len(next_frontier) >= size_limit:
                    yield next_frontier
                else:
                    stack.append(next_frontier)

if __name__ == '__main__':
    import sys
    import argparse
    import treeswift
    parser = argparse.ArgumentParser(description="Enumerate all subgraphs of a given size.")
    parser.add_argument('tree', help="Newick tree")
    parser.add_argument('--size', type=int, default=2, help="Size of subgraphs to enumerate")
    args = parser.parse_args()
    t = treeswift.read_tree_newick(args.tree)
    for subgraph in enumerate_subgraphs(t, args.size):
        print(subgraph)