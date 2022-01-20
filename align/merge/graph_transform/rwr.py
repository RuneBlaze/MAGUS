def rwrTransform(graph):
    from align.merge.magus_night import MagusNight
    return MagusNight.apply_transformation(graph.context, 'rwr')