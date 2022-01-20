using Base.Filesystem
function apply_transformation(c, trans)
    if trans == "rwr" 
        outf = joinpath(dirname(c.graph.graphPath), "graph_rwr.txt")
        if isfile(outf)
            return outf
        end
        g = AlnGraph(AlnContext(c.workingDir, c.subalignmentPaths); ghostrun = true)
        labels, adj = read_graph(c.graph.graphPath)
        newgraph = rwr_normalize_graph(labels, adj, g)
        open(outf, "w+") do f
            dump_graph_to_file(f, newgraph)
        end
        return outf
    else
        error("Unknown transformation: " * trans)
    end
end