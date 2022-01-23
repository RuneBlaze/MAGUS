using Base.Filesystem
include("RWR.jl")
include("MCL.jl")
function apply_transformation(c, trans)
    outf = joinpath(dirname(c.graph.graphPath), "graph_$(trans).txt")
    configs = c.configs
    if isfile(outf)
        return outf
    end
    if trans == "rwr"
        g = AlnGraph(AlnContext(c.workingDir, c.subalignmentPaths); ghostrun = true)
        labels, adj = read_graph(c.graph.graphPath)
        newgraph = rwr_normalize_graph(labels, adj, g)
    elseif trans == "mcl_filter"
        configs.log("mcl_filter will be run with I = $(configs.mclInflationFactor)")
        g = AlnGraph(AlnContext(c.workingDir, c.subalignmentPaths); ghostrun = true)
        labels, adj = read_graph(c.graph.graphPath)
        newgraph = mcl_filter(adj, c; inflation = configs.mclInflationFactor)
    elseif trans == "mcl_reweight"
        configs.log("mcl_reweight will be run with I = $(configs.mclInflationFactor) and its = $(configs.mclReweightIts)")
        newgraph = mcl_weight_replacement(c; inflation = configs.mclInflationFactor, iterations = configs.mclReweightIts)
    else
        error("Unknown transformation: " * trans)
    end
    open(outf, "w+") do f
        dump_graph_to_file(f, newgraph)
    end
    return outf
end
