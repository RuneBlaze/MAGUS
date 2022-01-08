using DataStructures
using SparseArrays
import Base.IdSet
using Base.Iterators

const TDict = OrderedDict

mutable struct AlnGraph # analog to AlignmentGraph, in Julia
    size :: Int
    matrix :: Vector{TDict{Int, Int}}
    clusters :: Vector{Vector{Int}}
    mat_subposmap :: Vector{Tuple{Int, Int}}
    subaln_lengths :: Vector{Int}
    subset_matrix_ind :: Vector{Int}

    workingdir :: String
    graphpath :: String
    clusterpath :: String
    tracepath :: String
end


# We really should not be using this struct, but we should keep things in lists/dicts
# I really think that using IdSet hurts performance here
mutable struct Node
    children :: IdSet{Node}
    parent :: Union{Nothing, Node}
    label :: Union{Nothing, Int}
    alive :: Bool
    num_elements :: Int
    rows :: BitSet
    outedges :: IdSet{Node}
    inedges :: IdSet{Node}
end

struct FlatCluster
    nodes :: Vector{Int}
    rows :: BitSet
    positions :: Vector{Tuple{Int, Int}}
end

struct ClusteringConfig
    prevent_vertical_violations :: Bool
    ensure_order :: Bool
    zero_weight :: Bool
    ClusteringConfig(pvv = true, eo = true, zw = false) = new(pvv, eo, zw)
end

# ClusteringConfig() = ClusteringConfig(true, true, true)
# ClusteringConfig(a = true, b = true, c = true) = ClusteringConfig(a, b, c)

function convert_to_flatclusters(clusters, alngraph :: AlnGraph)
    res = FlatCluster[]
    for cluster = clusters
        rows = BitSet()
        for c = cluster
            r = get_node_row(c, alngraph)
            @assert r ∉ rows
            push!(rows, r)
        end

        positions = Tuple{Int, Int}[]
        for c = cluster
            push!(positions, get_node_pos(c, alngraph))
        end
        push!(res, FlatCluster(cluster, rows, positions))
    end
    res
end

function partial_lessthan(lhs :: FlatCluster, rhs :: FlatCluster)
    for i = lhs.rows ∩ rhs.rows
        lhs_columns = [b for (a, b) = lhs.positions if a == i]
        rhs_columns = [b for (a, b) = rhs.positions if a == i]
        if maximum(lhs_columns) >= minimum(rhs_columns)
            return false
        end
    end
    return true
end

function partial_wellordered(lhs :: FlatCluster, rhs :: FlatCluster)
    if !isempty(lhs.rows ∩ rhs.rows)
        if partial_lessthan(lhs, rhs)
            return -1
        elseif partial_lessthan(rhs, lhs)
            return 1
        else
            error("Clusters are not well-ordered")
        end
    else
        return 0
    end
end

function check_flatclusters_validity(flatclusters)
    d = Dict(-1 => 0, 0 => 0, 1 => 0)
    n = length(flatclusters)
    for i = 1:n
        for j = i+1:n
            lhs = flatclusters[i]
            rhs = flatclusters[j]
            d[partial_wellordered(lhs, rhs)] += 1
        end
    end
    return d
end

function naive_newick(n :: Node)
    if !isnothing(n.label)
        return "$(n.label)"
    else
        mid = join([naive_newick(i) for i in n.children], ",")
        return "($mid)"
    end
end

Base.show(io::IO, x::Node) = print(io, naive_newick(x))

function Base.isless(lhs :: Node, rhs :: Node)
    if !isnothing(lhs.label) && !isnothing(rhs.label)
        return Base.isless(lhs.label, rhs.label)
    end
    return Base.isless(objectid(lhs), objectid(rhs))
end

function Node(a, b, c, d, e, f)
    return Node(a, b, c, d, e, f, IdSet(), IdSet())
end

function Node()
    return Node(Set(), nothing, nothing, true, 0, BitSet())
end

function Node(i :: Int)
    return Node(Set(), nothing, i, true, 1, BitSet())
end

@inline get_node_pos(i :: Int, g :: AlnGraph) = g.mat_subposmap[i + 1]
@inline get_node_row(i :: Int, g :: AlnGraph) = get_node_pos(i, g)[1]

function Node(i :: Int, g :: AlnGraph)
    return Node(Set(), nothing, i, true, 1, BitSet([get_node_row(i, g)]))
end

function connect!(u :: Node, v :: Node)
    push!(u.outedges, v)
    push!(v.inedges, u)
end

function disconnect!(u :: Node, v :: Node)
    delete!(u.outedges, v)
    delete!(v.inedges, u)
end

function exists_selfloop(u :: Node, v :: Node)
    return u ∈ v.outedges && v ∈ u.outedges
end

function contract!(u :: Node, v :: Node, mid :: Node)
    # mid is the parent node containing both u and v, we are to replace u and v with mid everywhere
    @assert !exists_selfloop(u, v)
    mid.outedges = u.outedges ∪ v.outedges
    mid.inedges = u.inedges ∪ v.inedges
    for e = mid.inedges
        delete!(e.outedges, u)
        delete!(e.outedges, v)
        push!(e.outedges, mid)
    end
    for e = mid.outedges
        delete!(e.inedges, u)
        delete!(e.inedges, v)
        push!(e.inedges, mid)
    end
end

function expand!(u :: Node, v :: Node, mid :: Node) # inverse of contract!
    for e = mid.inedges
        delete!(e.outedges, mid)
        delete!(e.outedges, mid)
    end

    for e = mid.outedges
        delete!(e.inedges, mid)
        delete!(e.inedges, mid)
    end

    for e = u.outedges
        push!(e.inedges, u)
    end

    for e = u.inedges
        push!(e.outedges, u)
    end

    for e = v.outedges
        push!(e.inedges, v)
    end

    for e = v.inedges
        push!(e.outedges, v)
    end
end

function dfs_exists_partial_order(root :: Node, u :: Node, v :: Node)
    if exists_selfloop(u, v)
        return false
    end

    function reachable_from(s, t :: Node)
        stack = collect(s)
        visited = Set{Node}()
        while !isempty(stack)
            n = pop!(stack)
            for e = n.outedges
                if e ∉ visited
                    if e == t
                        return true
                    end
                    push!(visited, e)
                    push!(stack, e)
                end
            end
        end
        return false
    end

    loop_exists = reachable_from(u.outedges, v) || reachable_from(v.outedges, u)
    return !loop_exists
end

function exists_partial_order(root :: Node, u :: Node, v :: Node)
    # check that after joining u and v in the tree rooted at root,
    # the partial order defined by incidence is still a valid partial order by toposort
    if exists_selfloop(u, v)
        return false
    end
    mid = Node() # intermediate node mid, that represents the join of u and v
    contract!(u, v, mid)
    # mid.outedges = u.outedges ∪ v.outedges

    # we do not modify the tree, but instead in toposort "swap out" u and v to be mid
    visited = IdSet{Node}()
    indegree = Dict{Node, Int}()
    queue = Queue{Node}()

    @inline function push_node!(node :: Node)
        if node ∉ visited
            push!(visited, node)
            enqueue!(queue, node)
        end
    end

    for c = root.children
        if c == u || c == v
            continue
        end
        indegree[c] = length(c.inedges)
        if indegree[c] == 0
            push_node!(c)
        end
    end

    indegree[mid] = length(mid.inedges)
    if indegree[mid] == 0
        push_node!(mid)
    end

    while !isempty(queue)
        # @show length(visited)
        node = dequeue!(queue)
        for e = node.outedges
            indegree[e] = indegree[e] - 1
            if indegree[e] == 0
                push_node!(e)
            end
        end
    end

    count = length(visited)

    @assert count <= length(root.children) - 1

    expand!(u, v, mid)
    # @show length(visited), length(root.children) - 1

    # total visited should be length(root.children) - 1
    if count < length(root.children) - 1
        # @show "$u $v rejected!"
        return false
    end

    return true
end

function add_child!(node :: Node, child :: Node)
    push!(node.children, child)
    child.parent = node
    node.num_elements += child.num_elements
end

@inline function numchildren(node :: Node)
    return length(node.children)
end

function remove_child!(node :: Node, child :: Node)
    delete!(node.children, child)
    node.num_elements -= child.num_elements
end

function join_nodes!(lhs :: Node, rhs :: Node)
    @assert lhs.parent == rhs.parent
    @assert isempty(lhs.rows ∩ rhs.rows)
    lhs.alive = false
    rhs.alive = false
    orig_parent = lhs.parent
    mid = Node()
    mid.rows = lhs.rows ∪ rhs.rows
    add_child!(mid, lhs)
    add_child!(mid, rhs)
    remove_child!(orig_parent, lhs)
    remove_child!(orig_parent, rhs)
    add_child!(orig_parent, mid)
    contract!(lhs, rhs, mid)
    # for e = mid.parent.children
    #     for c = e.outedges
    #         @assert c.alive
    #     end

    #     for c = e.inedges
    #         @assert c.alive "$(c == lhs) $(c == rhs) $(c == mid)"
    #     end
    # end
    mid
end

function is_valid_join(u :: Node, v :: Node; config :: ClusteringConfig)
    # @show config
    # error("showed config")
    root = u.parent
    if config.prevent_vertical_violations && !isempty(u.rows ∩ v.rows)
        return false
    end
    if config.ensure_order
        # old = exists_partial_order(root, u, v)
        new = dfs_exists_partial_order(root, u, v)
        # @assert old == new
        if !new
            return false
        end
    end
    return true
end

function star_tree(labels :: Vector{Int}, graph :: AlnGraph)
    leaves = [Node(n, graph) for n in labels]
    sort!(leaves; by = x -> x.label)
    bound = 0

    lengths = Iterators.Stateful(graph.subaln_lengths)

    total_connected = 0
    for i = 1:(length(leaves) - 1)
        first_node = leaves[i]
        second_node = leaves[i + 1]
        first_num = first_node.label
        second_num = second_node.label

        # invariant: the first num < bound
        while !(first_num < bound)
            bound += popfirst!(lengths)
        end
        if first_num < bound && second_num < bound
            connect!(first_node, second_node)
            total_connected += 1
        end
    end
    @show total_connected, graph.subaln_lengths, length(leaves), maximum(labels)
    root = Node()
    for l = leaves
        add_child!(root, l)
    end
    root
end

@inline function isleaf(node)
    return isempty(node.children)
end

function flat_labels(node :: Node)
    res = Int[]
    stack = []
    push!(stack, node)

    while !isempty(stack)
        top = pop!(stack)
        if isleaf(top)
            push!(res, top.label)
        else
            for e = top.children
                push!(stack, e)
            end
        end
    end
    return res
end

function upgma(labels :: Vector{Int}, similarity_ :: Dict{Int, Dict{Int, Float64}}, graph :: AlnGraph; 
    order = Base.Order.Reverse, config :: ClusteringConfig = ClusteringConfig())
    tree = star_tree(labels, graph)
    # pq = PriorityQueue{Tuple{Node, Node}, Float64}(order)
    pqs = Dict{Node, PriorityQueue{Node, Float64}}()
    for i=tree.children
        for j=tree.children
            if i < j && haskey(similarity_, i.label) && haskey(similarity_[i.label], j.label)
                # enqueue!(pq, (i, j), similarity_[i.label][j.label])
                
                enqueue!(get!(pqs, i, PriorityQueue(order)), j => similarity_[i.label][j.label])
            end
        end
    end
    return upgma_step2(tree, pqs; config = config)
end

function upgma_step2(tree :: Node, pqs :: Dict{Node, PriorityQueue{Node, Float64}}; config :: ClusteringConfig)
    @inline mkpair(a, b) = a < b ? (a, b) : (b, a)

    zero_weight_mode = config.zero_weight

    N = Threads.nthreads()

    while numchildren(tree) > 2
        thread_results = [Tuple{Float64, Node, Node}[] for _ = 1:N]
        dead_ls = [Node[] for _ = 1:N]
        Threads.@threads for (l, pq) = pqs
            tid = Threads.thread_id()
            if !l.alive
                push!(dead_ls[tid], l)
                continue
            end
            while !isempty(pq)
                r = peek(pq)
                v = pq[r]
                
                if !r.alive || !(is_valid_join(l, r; config))
                    dequeue!(pq)
                    continue
                end

                push!(thread_results[tid], (v, l, r))
                break
            end
        end

        if all(isempty, thread_results)
            break
        end

        (_, l, r) = maximum((v, l, r) for e = thread_results for (v, l, r) = e)
        dequeue!(pqs[l])
        dead_keys = [l for e = dead_ls for l = e]
        for e = dead_keys
            delete!(pqs, e)
        end

        mid = join_nodes!(l, r)

        for c = tree.children
            if c == mid || l == c || r == c
                continue
            end
            p1 = mkpair(l, c)
            p2 = mkpair(r, c)
            @inline in_pq(pos) = haskey(pqs, pos[1]) && haskey(pqs[pos[1]], pos[2])
            @inline get_pq(pos) = pqs[pos[1]][pos[2]]
            if in_pq(p1) && in_pq(p2)
                ud = (get_pq(p1) * l.num_elements + get_pq(p2) * r.num_elements)/(l.num_elements + r.num_elements)
            elseif in_pq(p1)
                if !zero_weight_mode
                    ud = get_pq(p1)
                else
                    ud = (get_pq(p1) * l.num_elements)/(l.num_elements + r.num_elements)
                end
            elseif in_pq(p2)
                if !zero_weight_mode
                    ud = get_pq(p2)
                else
                    ud = (get_pq(p2) * r.num_elements)/(l.num_elements + r.num_elements)
                end
            else
                continue
            end

            lhs, rhs = mkpair(mid, c)
            enqueue!(get!(pqs, lhs, PriorityQueue(Base.Order.Reverse)), rhs => ud)
            
            # enqueue!(pq, mkpair(mid, c), ud)
        end
    end
    return tree
end

function read_graph(io)
    @inline reorder(x, y) = x < y ? (x, y) : (y, x)
    graph = Dict{Int, Dict{Int, Float64}}()
    labels = Set{Int}()
    for line in eachline(io)
        _1, _2, _3 = [n for n = split(strip(line))]
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

function breakup_clusters(root :: Node, max_size :: Int)
    res = []
    stack = []
    push!(stack, root)
    while !isempty(stack)
        top = pop!(stack)
        if top.num_elements <= max_size
            push!(res, top)
            continue
        end
        for e = top.children
            push!(stack, e)
        end
    end
    return res
end

function upgma_naive_clustering(labels :: Vector{Int}, graph :: Dict{Int, Dict{Int, Float64}}, alngraph :: AlnGraph;
    config :: ClusteringConfig = ClusteringConfig())
    tree = upgma(labels, graph, alngraph; config = config)
    [flat_labels(c) for c = tree.children]
end

function find_clusters(c; config = ClusteringConfig())
    # c, alignment context
    g = AlnGraph(AlnContext(c.workingDir, c.subalignmentPaths); ghostrun = true)
    labels, adj = read_graph(c.graph.graphPath)
    upgma_naive_clustering(labels, adj, g; config = config)
end

# labels, G = read_graph("/home/lbq/research/datasets/UnalignFragTree/high_frag/1000M1/R0/unaligned_all.magus_take2_false/graph/graph.txt")
# open("res.txt", "w+") do f
#     r = upgma_naive_clustering(labels, G; size_limit = 25)
#     println(f, r)
# end

# using BenchmarkTools
# upgma_naive_clustering(
#     [1,2,3,4,5],
#     Dict(1 => Dict(2 => 17.0, 3 => 21.0, 4 => 31.0, 5 => 23.0),
#     2 => Dict(3 => 30.0, 4 => 34.0, 5 => 21.0),
#     3 => Dict(4 => 28.0, 5 => 39.0),
#     4 => Dict(5 => 43.0)); size_limit = 3
# )