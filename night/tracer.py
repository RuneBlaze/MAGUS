import sys
from os.path import dirname, realpath, join

sys.path.append(dirname(dirname(realpath(__file__))))

from configuration import Configs
from align.merge.graph_trace.min_clusters import minClustersSearch
from align.merge.graph_cluster.clean_clusters import purgeClusterViolations, purgeDuplicateClusters
from align.merge.alignment_graph import AlignmentGraph
from align.alignment_context import AlignmentContext
from glob import glob

Configs.workingDir = "/Users/lbq/Downloads/sandia_data/magus10krun_norecurse"

subalnPaths = glob(join(Configs.workingDir,"subalignments","*"))
subalnPaths.sort(key = lambda x: int(x.split("_")[-1].split(".")[0]))

from timeit import timeit
def critical():
    g = AlignmentGraph(
    AlignmentContext(
        workingDir=Configs.workingDir,
        subalignmentPaths=subalnPaths,))
    g.initializeMatrix()
    g.readGraphFromFile(g.graphPath)
    g.readClustersFromFile(g.clusterPath)
    purgeDuplicateClusters(g)
    purgeClusterViolations(g)
    minClustersSearch(g)
print(f"{timeit(lambda: critical(), number = 5) / 5=} seconds")