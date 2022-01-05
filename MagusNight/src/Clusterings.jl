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
    outedges :: IdSet{Node}
    inedges :: IdSet{Node}
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

function is_valid_join(u :: Node, v :: Node)
    root = u.parent
    if !isempty(u.rows ∩ v.rows)
        return false
    end
    if !exists_partial_order(root, u, v)
        return false
    end
    # @show "valid!"
    return true
end

function star_tree(labels :: Vector{Int}, graph :: AlnGraph)
    leaves = [Node(n, graph) for n in labels]
    # @show sort(collect(Set([first(n.rows) for n in leaves])))
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

function upgma(labels :: Vector{Int}, similarity_ :: Dict{Int, Dict{Int, Float64}}, graph :: AlnGraph; order = Base.Order.Reverse)
    tree = star_tree(labels, graph)
    pq = PriorityQueue{Tuple{Node, Node}, Float64}(order)
    for i=tree.children
        for j=tree.children
            if i < j && haskey(similarity_, i.label) && haskey(similarity_[i.label], j.label)
                enqueue!(pq, (i, j), similarity_[i.label][j.label])
            end
        end
    end
    return upgma_step2(tree, pq)
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

function upgma_step2(tree :: Node, pq :: PriorityQueue{Tuple{Node, Node}, Float64})
    @inline mkpair(a, b) = a < b ? (a, b) : (b, a)

    while numchildren(tree) > 2 &&! isempty(pq)
        l, r = dequeue!(pq)
        if !(l.alive && r.alive)
            continue
        end

        if !(is_valid_join(l, r))
            continue
        end
        # @show length(pq)
        mid = join_nodes!(l, r)

        for c = tree.children
            if c == mid || l == c || r == c
                continue
            end
            # this should be UPGMA*
            # if is_valid_join(mid, c)
            p1 = mkpair(l, c)
            p2 = mkpair(r, c)
            if haskey(pq, p1) && haskey(pq, p2)
                ud = (pq[p1] * l.num_elements + pq[p2] * r.num_elements)/(l.num_elements + r.num_elements)
            elseif haskey(pq, p1)
                ud = pq[p1]
            elseif haskey(pq, p2)
                ud = pq[p2]
            else
                continue
            end
            
            enqueue!(pq, mkpair(mid, c), ud)
            # end
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

function upgma_naive_clustering(labels :: Vector{Int}, graph :: Dict{Int, Dict{Int, Float64}}, alngraph :: AlnGraph)
    tree = upgma(labels, graph, alngraph)
    [flat_labels(c) for c = tree.children]
end

function find_clusters(c)
    # c, alignment context
    g = AlnGraph(AlnContext(c.workingDir, c.subalignmentPaths))
    # K = length(c.subalignmentPaths) # number of constraint alignments, upper bound of sizes
    labels, adj = read_graph(c.graph.graphPath)
    upgma_naive_clustering(labels, adj, g)
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