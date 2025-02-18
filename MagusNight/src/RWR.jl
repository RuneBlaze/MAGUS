# random walk with restart
using SparseArrays
using LinearAlgebra

function rwr_matrix_for(cc_::Set{Int}, digraph::Dict{Int,Dict{Int,Float64}})
    cc = collect(cc_)
    N = length(cc)
    ele2node = Vector{Int}(undef, length(cc))
    for i = 1:N
        ele2node[i] = cc[i]
    end
    node2ele = Dict{Int,Int}()
    for i = 1:N
        node2ele[ele2node[i]] = i
    end

    I = Int[]
    J = Int[]
    V = Float32[]
    cnt = length(digraph)
    for (u, adj) in digraph
        cnt -= 1
        if u ∉ cc_
            continue
        end
        for (v, w) in adj
            if v ∉ cc_
                continue
            end
            push!(I, node2ele[u])
            push!(J, node2ele[v])
            push!(V, w)
            push!(I, node2ele[v])
            push!(J, node2ele[u])
            push!(V, w)
        end
    end
    A = sparse(I, J, V, N, N)
    D = Diagonal(sum(A; dims = 1)[1, :])^(-0.5)
    A = D * A * D
    return (A = A, forward = node2ele, backward = ele2node)
end

function rwr(s::Int, A; delta = 0.5, max_iters = 10000)
    N = size(A)[1]
    e = zeros(Float32, N)
    e[s] = 1
    r = copy(e)
    for i = 1:max_iters
        rn = delta * A * r + (1 - delta) * e
        if norm(rn - r) < 1e-5
            break
        end
        r = rn
    end
    return r
end

function rwr_normalize_graph(
    labels::Vector{Int},
    digraph::Dict{Int,Dict{Int,Float64}},
    alngraph::AlnGraph,
)
    newgraph = Dict{Int,Dict{Int,Float64}}()
    ccs = connected_components(labels, digraph, alngraph)
    for (cc_ix, cc_) in enumerate(ccs)
        println("dealing with new connected component $cc_ix out of $(length(ccs))")
        cc = cc_.nodes
        cc_iter = collect(cc)
        sort!(cc_iter)
        rwr_setup = rwr_matrix_for(cc, digraph)
        println("obtained setup for rwr")
        N = length(cc)
        A = rwr_setup.A
        # sA = Symmetric(A)
        R = zeros(Float32, N, N)
        forward = rwr_setup.forward
        backward = rwr_setup.backward
        for (cnt, i) in enumerate(cc_iter)
            if cnt % 2000 == 0
                println("Progress: $cnt out of $N, roughly $(round(cnt / N * 100))%")
            end
            start_node = forward[i]
            R[start_node, :] = rwr(start_node, A)
        end
        R += transpose(R)
        max_r = maximum(R)
        partial_graph = Dict{Int,BinaryMaxHeap{Tuple{Float32,Int}}}()
        for u in cc_iter
            partial_graph[u] = BinaryMaxHeap{Tuple{Float32,Int}}()
            newgraph[u] = Dict{Int,Float64}()
        end
        for i = 1:N
            u = cc_iter[i]
            for j = i+1:N
                v = cc_iter[j]
                push!(partial_graph[u], (R[forward[u], forward[v]], v))
                push!(partial_graph[v], (R[forward[v], forward[u]], u))
            end
        end

        max_w = 0
        degs = counter(Int)
        for (u, adj) in digraph
            for (v, w) in adj
                degs[v] += 1
                degs[u] += 1
                max_w = max(max_w, w)
            end
        end

        @inline reorder(x, y) = x < y ? (x, y) : (y, x)

        for u in cc_iter
            cnt = 0
            while cnt < degs[u] * 2 && !isempty(partial_graph[u])
                (w, v) = pop!(partial_graph[u])
                x, y = reorder(u, v)
                newgraph[x][y] = w / max_r * max_w
                cnt += 1
            end
        end
    end
    return newgraph
end

function dump_graph_to_file(io, dict::Dict{Int,Dict{Int,Float64}})
    for (u, adj) in dict
        for (v, w) in adj
            if w == round(w)
                println(io, "$u $v $(Int(w))")
            else
                println(io, "$u $v $w")
            end
        end
    end
end

function elementary_digraph_stats(labels, digraph)
    num_nodes = length(labels)
    num_edges = 0
    for (k, v) in digraph
        num_edges += length(v)
    end
    return (num_nodes, num_edges)
end
