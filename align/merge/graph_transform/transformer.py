from configuration import Configs
import time
from align.merge.graph_transform.rwr import rwrTransform
from align.merge.magus_night import MagusNight

def transformGraph(graph):
    time1 = time.time()
    method = Configs.graphTransformMethod
    if not method:
        Configs.log("No alignment graph transformation requested..")
        return
    else:
        Configs.log(f"Applying {method} transformation...")
        ngraph_path = MagusNight.apply_transformation(graph.context, method)
        graph.graphPath = ngraph_path
    # elif method == "rwr":
    #     Configs.log("Applying RWR transformation..")
    #     ngraph_path = rwrTransform(graph)
    #     graph.graphPath = ngraph_path
    # else:
    #     raise Exception("Unknown graph transformation method: " + method)
    time2 = time.time()
    Configs.log("Transformed the graph in {} sec..".format(time2-time1))