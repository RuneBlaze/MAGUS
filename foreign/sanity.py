# sanity check of the Kruskal_MST thing
from foreign.Kruskal_MST import build_groups_MST

# tree :: DendropyTree
# subsets :: [[String]]
def build_MST(tree, subsets):
    grouping = {}
    for i, taxa in enumerate(subsets):
        for taxon in taxa:
            grouping[taxon] = i + 1
    subsets_tree = build_groups_MST(tree, grouping)
    return subsets_tree