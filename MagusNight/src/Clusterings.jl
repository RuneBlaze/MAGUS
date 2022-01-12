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

mutable struct Node
    children :: IdSet{Node}
    parent :: Union{Nothing, Node}
    label :: Union{Nothing, Int}
    alive :: Bool
    num_elements :: Int
    rows :: BitSet
    columns :: Vector{Int}
    outedges :: IdSet{Node}
    inedges :: IdSet{Node}
end

function dominated(a :: Node, b :: Node)
    # checks if a is dominated by b
    return all(a.columns[i] < b.columns[i] for i = a.rows ∩ b.rows)
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

    ClusteringConfig(a = true, b = true, c = false) = new(a, b, c)
end

function connected_components(labels :: Vector{Int}, digraph :: Dict{Int, Dict{Int, Float64}}, g :: AlnGraph)
    graph = Dict{Int, Vector{Int}}()
    for i = labels
        graph[i] = Int[]
    end

    for (u, dict) = digraph
        for (v, _) = dict
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
            push!(cc, (nodes = visited, rows = BitSet(get_node_row(i, g) for i = visited)))
        end
        setdiff!(unvisited, visited)
    end
    return cc
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
    return (d, n, extrema(length(i.nodes) for i = flatclusters))
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

function Node(a, b, c, d, e, f, g)
    return Node(a, b, c, d, e, f, g, IdSet(), IdSet())
end

function Node()
    return Node(Set(), nothing, nothing, true, 0, BitSet(), [])
end

function node_with_similar(n :: Node)
    return Node(Set(), nothing, nothing, true, 0, BitSet(), zeros(Int, length(n.columns)))
end

function Node(i :: Int)
    return Node(Set(), nothing, i, true, 1, BitSet(), [])
end

@inline get_node_pos(i :: Int, g :: AlnGraph) = g.mat_subposmap[i + 1]
@inline get_node_row(i :: Int, g :: AlnGraph) = get_node_pos(i, g)[1]

function initial_columns_for_node(i :: Int, g :: AlnGraph)
    pos = get_node_pos(i, g)
    k = length(g.subaln_lengths)
    r = zeros(Int, k)
    r[pos[1]] = pos[2]
    return r
end

function Node(i :: Int, g :: AlnGraph)
    return Node(Set(), nothing, i, true, 1, BitSet([get_node_row(i, g)]), initial_columns_for_node(i, g))
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

function dfs_exists_partial_order(root :: Node, u :: Node, v :: Node)
    if exists_selfloop(u, v)
        return false
    end
    loop_exists = reachable_from(u.outedges, v) || reachable_from(v.outedges, u)
    return !loop_exists
end

function constant_exists_partial_order(root :: Node, u :: Node, v :: Node)
    mid = node_with_similar(u)
    mid.outedges = u.outedges ∪ v.outedges
    mid.inedges = u.inedges ∪ v.inedges
    mid.rows = u.rows ∪ v.rows
    for i = u.rows
        mid.columns[i] = u.columns[i]
    end
    for i = v.rows
        mid.columns[i] = v.columns[i]
    end
    r = all(dominated(mid, o) for o = mid.outedges) && all(dominated(i, mid) for i = mid.inedges)
    return r
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
    mid = node_with_similar(lhs)
    mid.rows = lhs.rows ∪ rhs.rows
    for i = lhs.rows
        mid.columns[i] = lhs.columns[i]
    end
    for i = rhs.rows
        mid.columns[i] = rhs.columns[i]
    end
    add_child!(mid, lhs)
    add_child!(mid, rhs)
    remove_child!(orig_parent, lhs)
    remove_child!(orig_parent, rhs)
    add_child!(orig_parent, mid)
    contract!(lhs, rhs, mid)
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
        b = dfs_exists_partial_order(root, u, v)
        if !b
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

function disjoint_set_partial_order_exists(ds :: IntDisjointSets{Int64}, outedges :: Vector{Vector{Int}}, u :: Int, v :: Int)
    # assumption: u and v are current clusters -- they are roots of disjoint sets
    @inline redirect(x :: Int) = find_root(ds, x)
    if any(redirect(o) == v for o = outedges[u]) || any(redirect(o) == u for o = outedges[v])
        return false
    end
    function dfs(u :: Int, v :: Int)
        stack = Int[]
        visited = Set{Int}()
        for e = Set(redirect(x) for x = outedges[u])
            push!(stack, e)
        end
        while !isempty(stack)
            n = pop!(stack)
            for e_ = outedges[n]
                e = redirect(e_)
                if e ∉ visited
                    if e == v
                        return true
                    end
                    push!(stack, e)
                    push!(visited, e)
                end
            end
        end
        return false
    end
    if dfs(u, v) || dfs(v, u)
        return false
    end
    return true
end

@inline mkpair(a, b) = a < b ? (a, b) : (b, a)
# This time let's just do it right
function fast_upgma(labels :: Vector{Int}, similarity_ :: Dict{Int, Dict{Int, Float64}}, graph :: AlnGraph; 
    config :: ClusteringConfig = ClusteringConfig())
    sort!(labels)
    clusters = IntDisjointSets(length(labels))
    # we C++ now
    pq = PriorityQueue{Tuple{Int, Int}, Float64}(Base.Order.Reverse)
    rows = Vector{BitSet}(undef, length(labels))
    clustersizes = ones(Int, length(labels))
    node2initialcluster = Dict{Int, Int}()
    for (i, l) = enumerate(labels)
        node2initialcluster[l] = i
        rows[i] = BitSet([get_node_row(l, graph)])
    end
    weight_map = Dict{Int, Dict{Int, Float64}}() # similarity between two clusters
    for (u, map) = similarity_
        for (v, value) = map
            if u == v
                continue
            end
            lhs = node2initialcluster[u]
            rhs = node2initialcluster[v]
            lhs, rhs = mkpair(lhs, rhs)
            if lhs ∉ keys(weight_map)
                weight_map[lhs] = Dict{Int, Float64}()
            end
            if rhs ∉ keys(weight_map)
                weight_map[rhs] = Dict{Int, Float64}()
            end
            weight_map[lhs][rhs] = value
            weight_map[rhs][lhs] = value
            enqueue!(pq, (lhs, rhs), value)
        end
    end

    order_outedges = Vector{Int}[]
    # order_inedges = Vector{Int}[]
    for _ = 1:length(labels)
        push!(order_outedges, Int[])
        # push!(order_inedges, Int[])
    end
    bound = 0
    lengths = Iterators.Stateful(graph.subaln_lengths)
    total_connected = 0
    for i = 1:(length(labels)-1)
        first_num = labels[i]
        second_num = labels[i + 1]

        # invariant: the first num < bound
        while !(first_num < bound)
            bound += popfirst!(lengths)
        end
        if first_num < bound && second_num < bound
            push!(order_outedges[node2initialcluster[first_num]], node2initialcluster[second_num])
            total_connected += 1
        end
    end

    absorbed = Set{Int}()
    invalidated = Set{Tuple{Int, Int}}()
    while !isempty(pq)
        _, v = peek(pq)
        l, r = dequeue!(pq)
        l, r = mkpair(l, r)
        if l ∈ absorbed || r ∈ absorbed
            continue
        end

        if (l, r) ∈ invalidated || v != weight_map[l][r]
            continue
        end

        if !isempty(rows[l] ∩ rows[r]) || !disjoint_set_partial_order_exists(clusters, order_outedges, l, r)
            delete!(weight_map[l], r)
            delete!(weight_map[r], l)
            push!(invalidated, (l, r))
            continue
        end

        # now we merge l and r
        n = root_union!(clusters, l, r)
        m = l == n ? r : l # m is the cluster being merged
        push!(absorbed, m)

        # we update the order graph. This is not the most efficient way to do things
        # but we don't care for now. This obviously has more allocations
        # than necessary
        order_outedges[n] = order_outedges[n] ∪ order_outedges[m]
        # let's try a naive disjoint set thing

        # we update the weights. Remember, we are doing UPGMA
        for c = keys(weight_map[l]) ∪ keys(weight_map[r])
            # our job is to assign weight_map[c][n].
            # we should be able to do this in-place
            if haskey(weight_map[l], c) && haskey(weight_map[r], c)
                weight_map[n][c] = (weight_map[l][c] * clustersizes[l] + weight_map[r][c] * clustersizes[r])/(clustersizes[l] + clustersizes[r])
            elseif haskey(weight_map[l], c)
                weight_map[n][c] = weight_map[l][c]
            else
                weight_map[n][c] = weight_map[r][c]
            end
            weight_map[c][n] = weight_map[n][c]
        end

        # update the sizes
        clustersizes[n] = clustersizes[l] + clustersizes[r]
        rows[n] = rows[l] ∪ rows[r]
    end

    # naively export the lists
    final_clusters = Dict{Int, Vector{Int}}()
    for i = labels
        cid = find_root(clusters, node2initialcluster[i])
        if !haskey(final_clusters, cid)
            final_clusters[cid] = Vector{Int}()
        end
        push!(final_clusters[cid], i)
    end
    return values(final_clusters)
end

function upgma(labels :: Vector{Int}, similarity_ :: Dict{Int, Dict{Int, Float64}}, graph :: AlnGraph; 
    order = Base.Order.Reverse, config :: ClusteringConfig = ClusteringConfig())
    tree = star_tree(labels, graph)
    pq = PriorityQueue{Tuple{Node, Node}, Float64}(order)
    for i=tree.children
        for j=tree.children
            if i < j && haskey(similarity_, i.label) && haskey(similarity_[i.label], j.label)
                enqueue!(pq, (i, j), similarity_[i.label][j.label])
            end
        end
    end
    return upgma_step2(tree, pq; config = config)
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

function upgma_step2(tree :: Node, pq :: PriorityQueue{Tuple{Node, Node}, Float64}; config :: ClusteringConfig)

    zero_weight_mode = config.zero_weight

    while numchildren(tree) > 2 &&! isempty(pq)
        l, r = dequeue!(pq)
        if !(l.alive && r.alive)
            continue
        end

        if !(is_valid_join(l, r; config))
            continue
        end
        mid = join_nodes!(l, r)
        for c = tree.children
            if c == mid || l == c || r == c
                continue
            end
            p1 = mkpair(l, c)
            p2 = mkpair(r, c)
            if haskey(pq, p1) && haskey(pq, p2)
                ud = (pq[p1] * l.num_elements + pq[p2] * r.num_elements)/(l.num_elements + r.num_elements)
            elseif haskey(pq, p1)
                if !zero_weight_mode
                    ud = pq[p1]
                else
                    ud = (pq[p1] * l.num_elements)/(l.num_elements + r.num_elements)
                end
            elseif haskey(pq, p2)
                if !zero_weight_mode
                    ud = pq[p2]
                else
                    ud = (pq[p2] * r.num_elements)/(l.num_elements + r.num_elements)
                end
            else
                continue
            end
            
            enqueue!(pq, mkpair(mid, c), ud)
        end
    end
    return tree
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