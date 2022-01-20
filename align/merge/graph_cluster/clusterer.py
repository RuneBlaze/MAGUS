'''
Created on Aug 23, 2020

@author: Vlad
'''

import os
import time

from configuration import Configs

from align.merge.graph_cluster.mcl import runMclClustering
from align.merge.graph_cluster.mlr_mcl import runMlrMclClustering
from align.merge.graph_cluster.rg import rgClustering

'''
The alignment graph is clustered, the clusters are written out as an array of node arrays.
MCL is the main way to do this, but rg could be used if there are scalability issues.
'''

def clusterGraph(graph):
    time1 = time.time()
    Configs.log(f"Received alignment graph with path {graph.graphPath}")
    
    if os.path.exists(graph.clusterPath):
        Configs.log("Found existing cluster file {}".format(graph.clusterPath))
        graph.readClustersFromFile(graph.clusterPath)
        
    elif Configs.graphClusterMethod == "mcl":
        runMclClustering(graph)
    elif Configs.graphClusterMethod == "upgma":
        from align.merge.magus_night import MagusNight
        upgma_config = MagusNight.ClusteringConfig(True, Configs.upgmaKeepOrder, Configs.upgmaZeroWeight, Configs.exp)
        Configs.log(f"Running UPGMA clustering/tracing with configuration: {upgma_config}")
        r = MagusNight.find_clusters(graph.context, config = upgma_config)
        graph.clusters = [list(e) for e in r]
        graph.writeClustersToFile(graph.clusterPath)
        graph.writeClustersToFile(graph.tracePath)
    elif Configs.graphClusterMethod == "mlrmcl":
        runMlrMclClustering(graph)
        
    elif Configs.graphClusterMethod == "rg":
        rgClustering(graph)
        
    else:
        Configs.log("No alignment graph clustering requested..")
    
    time2 = time.time()  
    Configs.log("Clustered the graph in {} sec..".format(time2-time1))
