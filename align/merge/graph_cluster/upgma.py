from configuration import Configs
from align.merge.graph_trace.upgma_star import upgmaCluster

def upgmaClustering(graph, zeroWeight = False):
    Configs.log("Building a UPGMA graph clustering with missing data strategy: {}..".format(
        "zero weight" if zeroWeight else "existing weight"))
    
    k = len(graph.context.subalignments) or len(graph.context.subalignmentPaths)
    lowerBound = [graph.subsetMatrixIdx[i] for i in range(k)]
    upperBound = [graph.subsetMatrixIdx[i] + graph.subalignmentLengths[i] for i in range(k)] 
    graph.clusters =  upgmaCluster(graph, lowerBound, upperBound, zeroWeight)
    graph.writeClustersToFile(graph.clusterPath)