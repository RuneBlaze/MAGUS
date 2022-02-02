from configuration import Configs
import time
from align.merge.graph_transform.rwr import rwrTransform

def transformGraph(graph):
    time1 = time.time()
    method = Configs.graphTransformMethod
    if not method:
        Configs.log("No alignment graph transformation requested..")
        return
    else:
        from align.merge.magus_night import MagusNight
        Configs.log(f"Applying {method} transformation...")
        ngraph_path = MagusNight.apply_transformation(graph.context, method)
        graph.graphPath = ngraph_path
        graph.readGraphFromFile(ngraph_path)
    time2 = time.time()
    Configs.log("Transformed the graph in {} sec..".format(time2-time1))