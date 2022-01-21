function mcl_filter(digraph :: Dict{Int, Dict{Int, Float64}}, context; inflation = 1.4)
    config = context.configs
    mcl_path = config.mclPath
    gp = context.graph.graphPath
    op = "$gp.mcl_filter"
    run(`$mcl_path $gp --abc -I $inflation -o $op`)
    clusters = Vector{Vector{Int}}()
    for line in eachline(op)
        tokens = [parse(Int, token) for token in split(strip(line))]
        push!(clusters, tokens)
    end

    cid = Dict{Int, Int}()
    for (i, c) = enumerate(clusters)
        for n = c
            cid[n] = i
        end
    end

    newgraph = Dict{Int, Dict{Int, Float64}}()

    for (u, adj) = digraph
        for (v, w) = adj
            if cid[u] == cid[v]
                if !haskey(newgraph, u)
                    newgraph[u] = Dict{Int, Float64}()
                end
                newgraph[u][v] = w
            end
        end
    end
    return newgraph
end

function mcl_weight_replacement(context;
    inflation = 2, iterations = 2)
    it = iterations
    config = context.configs
    mcl_path = config.mclPath
    gp = context.graph.graphPath
    op = "$gp.mcl_filter"
    # @show `$mcl_path $gp --abc -I $inflation -o $op -dump-interval $(it):$(it+1) -dump linesite`
    run(`$mcl_path $gp --abc -I $inflation --d -dump-stem dump -o $op -dump-interval $(it):$(it+1) -dump linesite`)
    graph_file = "ite-$it.dump"
    labels, new_graph = read_graph(graph_file)
    return new_graph
end