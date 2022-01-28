using DataStructures
using SparseArrays
import Base.IdSet
using Base.Iterators

const TDict = OrderedDict

mutable struct AlnGraph # analog to AlignmentGraph, in Julia
    size::Int
    matrix::Vector{TDict{Int,Int}}
    clusters::Vector{Vector{Int}}
    mat_subposmap::Vector{Tuple{Int,Int}}
    subaln_lengths::Vector{Int}
    subset_matrix_ind::Vector{Int}

    workingdir::String
    graphpath::String
    clusterpath::String
    tracepath::String
end

struct ClusteringConfig
    prevent_vertical_violations::Bool
    ensure_order::Bool
    zero_weight::Bool
    clean_forbidden::Bool
    rg_fast :: Bool

    ClusteringConfig(a = true, b = true, c = false, d = false, e = false) = new(a, b, c, d, e)
end

function read_graph(io)
    @inline reorder(x, y) = x < y ? (x, y) : (y, x)
    graph = Dict{Int,Dict{Int,Float64}}()
    labels = Set{Int}()
    for line in eachline(io)
        _1, _2, _3 = [n for n in split(strip(line))]
        a_, b_, w = parse(Int, _1), parse(Int, _2), parse(Float64, _3)
        a, b = reorder(a_, b_)
        if !haskey(graph, a)
            graph[a] = Dict()
        end
        push!(labels, a)
        push!(labels, b)
        graph[a][b] = w
    end
    return collect(labels), graph
end

function symmetrize!(g)
    for (u, adj) = g
        for (v, w) = adj
            if !haskey(g, v)
                g[v] = Dict{Int, eltype(w)}()
            end
            g[v][u] = w
        end
    end
end

function desymmetrize(g)
    newgraph = Dict{Int,Dict{Int,Float64}}()
    for (u, adj) = g
        for (v, w) = adj
            if u > v
                continue
            end
            if !haskey(g, u)
                g[u] = Dict{Int, eltype(w)}()
            end
            g[u][v] = w
        end
    end
    return newgraph
end


include("ClusteringsOld.jl")

function connected_components(
    labels::Vector{Int},
    digraph::Dict{Int,Dict{Int,Float64}},
    g::AlnGraph,
)
    graph = Dict{Int,Vector{Int}}()
    for i in labels
        graph[i] = Int[]
    end

    for (u, dict) in digraph
        for (v, _) in dict
            if u != v
                push!(graph[u], v)
                push!(graph[v], u)
            end
        end
    end

    function dfs(node)
        visited = Set{Int}()
        stack = Vector{Int}()
        push!(stack, node)
        while !isempty(stack)
            n = pop!(stack)
            for v in graph[n]
                if v ∉ visited
                    push!(stack, v)
                    push!(visited, v)
                end
            end
        end
        visited
    end

    unvisited = Set{Int}(labels)
    cc = []
    while !isempty(unvisited)
        n = pop!(unvisited)
        visited = dfs(n)
        if !isempty(visited)
            push!(cc, (nodes = visited, rows = BitSet(get_node_row(i, g) for i in visited)))
        end
        setdiff!(unvisited, visited)
    end
    return cc
end

@inline get_node_pos(i::Int, g::AlnGraph) = g.mat_subposmap[i+1]
@inline get_node_row(i::Int, g::AlnGraph) = get_node_pos(i, g)[1]

function ds_dfs(outedges::Vector{Vector{Int}}, visited::BitVector, u::Int, v::Int)
    stack = Int[]
    for e in Set(x for x in outedges[u])
        push!(stack, e)
    end
    while !isempty(stack)
        n = pop!(stack)
        for e in outedges[n]
            if !visited[e]
                if e == v
                    return true
                end
                push!(stack, e)
                visited[e] = true
            end
        end
    end
    return false
end

function disjoint_set_partial_order_exists(
    visited::BitVector,
    outedges::Vector{Vector{Int}},
    u::Int,
    v::Int,
    mid::Int,
    mid_visited::BitVector,
)
    if v ∈ outedges[u] && u ∈ outedges[v]
        return false
    end

    @inline dfs(a, b, c) = ds_dfs(outedges, a, b, c)

    fill!(visited, false)
    if mid_visited[u] && mid_visited[v]
        if dfs(visited, u, v) || dfs(visited, v, u)
            return false
        end
    elseif mid_visited[u] # mid < u, v not
        if dfs(visited, v, u)
            return false
        end
    elseif mid_visited[v]
        if dfs(visited, u, v)
            return false
        end
    else
        if dfs(visited, u, v) || dfs(visited, v, u)
            return false
        end
    end

    return true
end

function order_clusters(absorbed::BitVector, outedges::Vector{Vector{Int}})
    N = length(outedges)
    indegs = zeros(Int, N)
    for (ix, i) in enumerate(outedges)
        if !absorbed[ix]
            for e in i
                indegs[e] += 1
            end
        end
    end

    stack = Int[]
    for (i, e) in enumerate(indegs)
        if !absorbed[i] && e == 0
            push!(stack, i)
        end
    end

    order = Int[]
    # naive toposort
    while !isempty(stack)
        t = pop!(stack)
        push!(order, t)
        for e in outedges[t]
            indegs[e] -= 1
            if indegs[e] == 0
                push!(stack, e)
            end
        end
    end
    return order
end

@inline mkpair(a, b) = a < b ? (a, b) : (b, a)
# This time let's just do it right
function fast_upgma(
    labels::Vector{Int},
    similarity_::Dict{Int,Dict{Int,Float64}},
    graph::AlnGraph;
    config::ClusteringConfig = ClusteringConfig(),
)
    sort!(labels)
    N = length(labels)
    clusters = IntDisjointSets(length(labels))
    # we C++ now
    # pq = BinaryMaxHeap{Tuple{Float64, Int, Int}}()
    pq = BinaryHeap{Tuple{Float64,Int,Int},DataStructures.FasterReverse}()
    rows = Vector{BitSet}(undef, length(labels))
    clustersizes = ones(Int, length(labels))
    node2initialcluster = Dict{Int,Int}()
    weight_map = Vector{Dict{Int,Float64}}(undef, length(labels))
    visited = BitVector(undef, length(labels))
    mid_visited = falses(length(labels))
    for (i, l) in enumerate(labels)
        node2initialcluster[l] = i
        rows[i] = BitSet([get_node_row(l, graph)])
        weight_map[i] = Dict{Int,Float64}()
    end

    for (u, map) in similarity_
        for (v, value) in map
            if u == v
                continue
            end
            lhs = node2initialcluster[u]
            rhs = node2initialcluster[v]
            lhs, rhs = mkpair(lhs, rhs)
            weight_map[lhs][rhs] = value
            weight_map[rhs][lhs] = value
            push!(pq, (value, lhs, rhs))
        end
    end

    order_outedges = Vector{Vector{Int}}(undef, length(labels))
    order_inedges = Vector{Vector{Int}}(undef, length(labels))
    for i = 1:length(labels)
        order_outedges[i] = Int[]
        order_inedges[i] = Int[]
    end
    bound = 0
    lengths = Iterators.Stateful(graph.subaln_lengths)
    total_connected = 0
    for i = 1:(length(labels)-1)
        first_num = labels[i]
        second_num = labels[i+1]

        # invariant: the first num < bound
        while !(first_num < bound)
            bound += popfirst!(lengths)
        end
        if first_num < bound && second_num < bound
            push!(
                order_outedges[node2initialcluster[first_num]],
                node2initialcluster[second_num],
            )
            push!(
                order_inedges[node2initialcluster[second_num]],
                node2initialcluster[first_num],
            )
            total_connected += 1
        end
    end

    absorbed = falses(N)
    invalidated = Set{Tuple{Int,Int}}()
    mid = -1

    @inline isforbidden(l, r) = !isempty(rows[l] ∩ rows[r]) || !disjoint_set_partial_order_exists(
            visited,
            order_outedges,
            l,
            r,
            mid,
            mid_visited,
        )
    @inline function forbid(l, r)
        delete!(weight_map[l], r)
        delete!(weight_map[r], l)
        push!(invalidated, (l, r))
    end
    while !isempty(pq)
        v, l, r = pop!(pq)
        l, r = mkpair(l, r)
        if absorbed[l] || absorbed[r]
            continue
        end

        if (l, r) ∈ invalidated || v != weight_map[l][r]
            continue
        end

        if isforbidden(l, r)
            forbid(l, r)
            continue
        end

        mid_visited = visited
        visited = BitVector(undef, N)
        # now we merge l and r
        n = root_union!(clusters, l, r)
        m = l == n ? r : l # m is the cluster being merged
        absorbed[m] = true

        if config.clean_forbidden
            # filter out edges that are no longer valid
            for c = keys(weight_map[l])
                _1, _2 = mkpair(l, c)
                if isforbidden(_1, _2)
                    forbid(_1, _2)
                end
            end
            for c = keys(weight_map[r])
                _1, _2 = mkpair(r, c)
                if isforbidden(_1, _2)
                    forbid(_1, _2)
                end
            end
        end

        # perform contraction
        union!(order_outedges[n], order_outedges[m])
        union!(order_inedges[n], order_inedges[m])
        for innode in order_inedges[m]
            for (i, e) in enumerate(order_outedges[innode])
                if e == m
                    order_outedges[innode][i] = n
                end
            end
        end
        for outnode in order_outedges[m]
            for (i, e) in enumerate(order_inedges[outnode])
                if e == m
                    order_inedges[outnode][i] = n
                end
            end
        end

        # we update the weights. Remember, we are doing UPGMA
        for c in keys(weight_map[l]) ∪ keys(weight_map[r])
            # our job is to assign weight_map[c][n].
            # we should be able to do this in-place
            if config.rg_fast
                weight_map[n][c] = get(weight_map[l], c, 0) + get(weight_map[r], c, 0)
            elseif config.zero_weight
                weight_map[n][c] =
                    (
                        get(weight_map[l], c, 0) * clustersizes[l] +
                        get(weight_map[r], c, 0) * clustersizes[r]
                    ) / (clustersizes[l] + clustersizes[r])
            else
                if haskey(weight_map[l], c) && haskey(weight_map[r], c)
                    weight_map[n][c] =
                        (
                            weight_map[l][c] * clustersizes[l] +
                            weight_map[r][c] * clustersizes[r]
                        ) / (clustersizes[l] + clustersizes[r])
                elseif haskey(weight_map[l], c)
                    weight_map[n][c] = weight_map[l][c]
                else
                    weight_map[n][c] = weight_map[r][c]
                end
            end
            weight_map[c][n] = weight_map[n][c]
            push!(pq, (weight_map[c][n], mkpair(c, n)...))
        end

        # update the sizes
        clustersizes[n] += clustersizes[m]
        mid = n
        union!(rows[n], rows[m])
    end

    cluster_ordering = order_clusters(absorbed, order_outedges)
    # naively export the lists
    final_clusters = Dict{Int,Vector{Int}}()
    for i in labels
        cid = find_root(clusters, node2initialcluster[i])
        if !haskey(final_clusters, cid)
            final_clusters[cid] = Vector{Int}()
        end
        push!(final_clusters[cid], i)
    end
    clusters = Vector{Vector{Int}}(undef, length(cluster_ordering))
    for (i, j) in enumerate(cluster_ordering)
        clusters[i] = final_clusters[j]
    end
    return clusters
end

function find_clusters(c; config = ClusteringConfig())
    # c, alignment context
    g = AlnGraph(AlnContext(c.workingDir, c.subalignmentPaths); ghostrun = true)
    labels, adj = read_graph(c.graph.graphPath)
    fast_upgma(labels, adj, g; config = config)
end
