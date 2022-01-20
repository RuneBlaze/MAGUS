from configuration import Configs
import time
from align.merge.graph_transform.rwr import rwrTransform

def transformGraph(graph):
    time1 = time.time()
    method = Configs.graphTransformMethod
    if not method:
        Configs.log("No alignment graph transformation requested..")
        return
    elif method == "rwr":
        Configs.log("Applying RWR transformation..")
        ngraph_path = rwrTransform(graph)
        graph.graphPath = ngraph_path
        # graph.readGraphFromFile(ngraph_path)
    else:
        raise("Unknown graph transformation method: " + method)
    time2 = time.time()
    Configs.log("Transformed the graph in {} sec..".format(time2-time1))