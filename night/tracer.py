import sys
from os.path import dirname, realpath, join

sys.path.append(dirname(dirname(realpath(__file__))))

from configuration import Configs
from align.merge.graph_trace.min_clusters import minClustersSearch
from align.merge.graph_cluster.rg import rgClustering
from align.merge.graph_cluster.mcl import runMclClustering
from align.merge.graph_cluster.clean_clusters import purgeClusterViolations, purgeDuplicateClusters
from align.merge.alignment_graph import AlignmentGraph
from align.alignment_context import AlignmentContext
from glob import glob
import julia

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
    # minClustersSearch(g)
    g.writeClustersToFile("scratch/clusters.python.txt")

from julia import Julia
jl = Julia()
jl.eval('push!(LOAD_PATH, "MagusNight/")')

from julia import MagusNight

def julia_critical():
    c = AlignmentContext(
        workingDir=Configs.workingDir,
        subalignmentPaths=subalnPaths,)
    g = AlignmentGraph(c)
    g.initializeMatrix()
    g.readGraphFromFile(g.graphPath)
    g.readClustersFromFile(g.clusterPath)
    results = MagusNight.find_trace(c)
    clusters = []
    for c in results.clusters:
        clusters.append(list(c - 1))
    g.clusters = clusters

    g.writeClustersToFile("scratch/clusters.night2.txt")

def julia_clustering():
    c = AlignmentContext(
        workingDir=Configs.workingDir,
        subalignmentPaths=subalnPaths,
        subalignments = subalnPaths)
    g = AlignmentGraph(c)
    c.graph = g
    g.context.subalignments = subalnPaths
    g.initializeMatrix()
    g.readGraphFromFile(g.graphPath)
    g.readClustersFromFile(g.clusterPath)
    results = MagusNight.find_clusters(c, config = MagusNight.ClusteringConfig(True, True))
    g.clusters = [list(e) for e in results]
    purgeDuplicateClusters(g)
    purgeClusterViolations(g)
    minClustersSearch(g)
    g.writeClustersToFile("scratch/ordered_clusters.txt")
    # print(results)

def python_rg_clustering():
    c = AlignmentContext(
        workingDir=Configs.workingDir,
        subalignmentPaths=subalnPaths,)
    g = AlignmentGraph(c)
    c.graph = g
    
    g.initializeMatrix()
    g.readGraphFromFile(g.graphPath)

    # g.readClustersFromFile(g.clusterPath)
    rgClustering(g, False)
    # runMclClustering(g)
    print(f"{g.clusters=}")
    tt = sum(len(c) for c in g.clusters)
    print(f"{tt=}")
    print(MagusNight.check_flatclusters_validity(MagusNight.convert_to_flatclusters(g.clusters, MagusNight.AlnGraph(c))))

# critical()
# julia_critical()
# julia_clustering()
python_rg_clustering()
# print(f"{timeit(lambda: critical(), number = 5) / 5=} seconds")