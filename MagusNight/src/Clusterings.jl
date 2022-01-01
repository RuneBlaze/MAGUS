using DataStructures
using SparseArrays
import Base.IdSet

mutable struct Node
    children :: IdSet{Node}
    parent :: Union{Nothing, Node}
    label :: Union{Nothing, Int}
    alive :: Bool
    num_elements :: Int
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

function Node()
    return Node(Set(), nothing, nothing, true, 0)
end

function Node(i :: Int)
    return Node(Set(), nothing, i, true, 1)
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
    lhs.alive = false
    rhs.alive = false
    orig_parent = lhs.parent
    mid = Node()
    add_child!(mid, lhs)
    add_child!(mid, rhs)
    remove_child!(orig_parent, lhs)
    remove_child!(orig_parent, rhs)
    add_child!(orig_parent, mid)
    mid
end

function star_tree(labels :: Vector{Int})
    leaves = [Node(n) for n in labels]
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

function upgma(labels :: Vector{Int}, similarity_ :: Dict{Int, Dict{Int, Float64}}; order = Base.Order.Reverse, size_limit = -1)
    tree = star_tree(labels)
    pq = PriorityQueue{Tuple{Node, Node}, Float64}(order)
    for i=tree.children
        for j=tree.children
            if i < j && haskey(similarity_, i.label) && haskey(similarity_[i.label], j.label)
                enqueue!(pq, (i, j), similarity_[i.label][j.label])
            end
        end
    end
    return upgma_step2(tree, pq; size_limit = size_limit)
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

function upgma_step2(tree :: Node, pq :: PriorityQueue{Tuple{Node, Node}, Float64}; size_limit :: Int)
    @inline mkpair(a, b) = a < b ? (a, b) : (b, a)

    while numchildren(tree) > 2 &&! isempty(pq)
        l, r = dequeue!(pq)
        if !(l.alive && r.alive)
            continue
        end
        mid = join_nodes!(l, r)

        for c = tree.children
            if c == mid || l == c || r == c
                continue
            end
            # this should be UPGMA*
            if size_limit < 0 || mid.num_elements + c.num_elements <= size_limit
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
            end
        end

        # @inline function remove_from_pq!(a, b)
        #     p = mkpair(a, b)
        #     if p in keys(pq)
        #         delete!(pq, p)
        #     end
        # end

        # for i = tree.children
        #     remove_from_pq!(l, i)
        #     remove_from_pq!(r, i)
        # end
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

function upgma_naive_clustering(labels :: Vector{Int}, graph :: Dict{Int, Dict{Int, Float64}}; size_limit :: Int = -1)
    tree = upgma(labels, graph; size_limit = size_limit)
    [flat_labels(c) for c = tree.children]
end

function find_clusters(c)
    # c, alignment context
    K = length(c.subalignmentPaths) # number of constraint alignments, upper bound of sizes
    labels, adj = read_graph(c.graph.graphPath)
    upgma_naive_clustering(labels, adj; size_limit = K)
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