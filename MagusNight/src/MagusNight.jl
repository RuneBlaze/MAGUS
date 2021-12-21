# reads a graph, and outputs a valid trace
module MagusNight
using SparseArrays
using Base.Filesystem
using DataStructures
using Debugger
# using Parsers

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

mutable struct AlnContext
    workingdir :: String
    subalnpaths :: Vector{String}
end

function read_seqlen_from_fasta(filename :: AbstractString)
    len = 0
    read_seq = false
    for line in eachline(filename)
        line = strip(line)
        if startswith(line, ">")
            if read_seq
                return len
            end
            read_seq = true
        else
            len += length(line)
        end
    end
    if read_seq
        return len
    else
        error("no sequence found in fasta file")
    end
end

function from_python_alngraph(pgraph)
    from = pgraph
    alnlengths = pgraph.subalignmentLengths

    size = sum(alnlengths)
    @show size, alnlengths
    subset_matrix_ind = zeros(Int, length(alnlengths))
    for k in 2:length(alnlengths)
        # some sort of prefix sum that I don't understand
        subset_matrix_ind[k] = subset_matrix_ind[k-1] + alnlengths[k-1]
    end
    
    mat_subposmap = Vector{Tuple{Int, Int}}(undef, size)
    i = 1
    for k in 1:length(alnlengths)
        for j in 1:alnlengths[k]
            mat_subposmap[i] = (k, j)
            i += 1
        end
    end

    pymatrix = from.matrix
    matrix = Vector{TDict{Int, Int}}(undef, size)
    for i = 1:size
        matrix[i] = TDict{Int, Int}()
    end

    for i in keys(pymatrix)
        for j in keys(pymatrix[i])
            matrix[i+1][j+1] = pymatrix[i+1][j+1]
            # s += pymatrix[i][j]
        end
    end

    clusters = Vector{Vector{Int}}()
    pyclusters = from.clusters
    for c in pyclusters
        push!(clusters, c .+ 1)
    end

    return AlnGraph(size, matrix, clusters, mat_subposmap,
        alnlengths, subset_matrix_ind, from.workingDir, from.graphPath, from.clusterPath, from.tracePath)
end

function AlnGraph(context :: AlnContext)
    wp = joinpath(context.workingdir, "graph")
    gp = joinpath(wp, "graph.txt")
    cp = joinpath(wp, "clusters.txt")
    tp = joinpath(wp, "trace.txt")

    for p in [gp, cp, tp]
        if !isfile(p)
            error("File not found: " * p)
        end
    end

    alnlengths = Int[]
    for (i, p) in enumerate(context.subalnpaths)
        push!(alnlengths, read_seqlen_from_fasta(p))
    end

    size = sum(alnlengths)
    # @show size, alnlengths
    subset_matrix_ind = zeros(Int, length(alnlengths))
    for k in 2:length(alnlengths)
        # some sort of prefix sum that I don't understand
        subset_matrix_ind[k] = subset_matrix_ind[k-1] + alnlengths[k-1]
    end
    
    mat_subposmap = Vector{Tuple{Int, Int}}(undef, size)
    i = 1
    for k in 1:length(alnlengths)
        for j in 1:alnlengths[k]
            mat_subposmap[i] = (k, j)
            i += 1
        end
    end
    matrix = Vector{TDict{Int, Int}}(undef, size)
    for i = 1:size
        matrix[i] = TDict{Int, Int}()
    end

    for line in eachline(gp)
        a, b, c = [parse(Int, token) for token in split(strip(line))]
        matrix[a+1][b+1] = c
        # s += c
    end

    clusters = Vector{Vector{Int}}()
    for line in eachline(cp)
        tokens = [parse(Int, token) + 1 for token in split(strip(line))]
        if length(tokens) > 1
            push!(clusters, tokens)
        end
    end

    # @show length(clusters), sum([sum(c) for c in clusters])

    return AlnGraph(size, matrix, clusters, mat_subposmap,
        alnlengths, subset_matrix_ind, wp, gp, cp, tp)
end

State = @NamedTuple begin
    heuristic :: Tuple{Float64, Int, Int, Int}
    numordered :: Int
    numleft :: Int
    pairsleft :: Int
    counter :: Int
    queue_idxs :: TDict{Int, Int}
    cluster_breaks :: TDict{Tuple{Int, Int}, Vector{Int}}
    maximal_cut :: TDict{Int, Int}
    newcluster_breaks :: Vector{Tuple{Int, Vector{Int}, Vector{Int}, Set{Int}}}
    safefrontier :: Bool
end

function Base.empty!(heap :: AbstractHeap)
    Base.empty!(heap.valtree)
end

function purge_duplicate_clusters(graph :: AlnGraph)
    unique_clusters = Set()
    new_clusters = []
    for cluster = graph.clusters
        sort!(cluster)
        cluster_tuple = (cluster...,)
        if cluster_tuple ∉ unique_clusters
            push!(unique_clusters, cluster_tuple)
            push!(new_clusters, cluster)
        end
    end
    graph.clusters = new_clusters
    println("Purged duplicate clusters. Found $(length(graph.clusters)) unique clusters..")
end

function purge_cluster_violations(graph :: AlnGraph)
    redundant_cols = TDict()
    redundant_rows = TDict()
    element_scores = TDict()

    for (a, cluster) = enumerate(graph.clusters)
        for b = cluster
            bsub, bpos = graph.mat_subposmap[b]
            redundant_cols[(a, bsub)] = get(redundant_cols, (a, bsub), Set()) ∪ Set([(a, b)])
            redundant_rows[b] = get(redundant_rows, b, Set()) ∪ Set([(a, b)])

            scoresum = 0
            for c = cluster
                csub, cpos = graph.mat_subposmap[c]
                if bsub != csub
                    scoresum += get(graph.matrix[b], c, 0)
                end
            end
            element_scores[(a, b)] = scoresum
        end
    end

    problem_cols = [(a, b) for (a, b) in keys(redundant_cols) if length(redundant_cols[(a, b)]) > 1]
    problem_rows = [a for a in keys(redundant_rows) if length(redundant_rows[a]) > 1]
    println("Found $(length(problem_rows)) row violations and $(length(problem_cols)) column violations")

    sorted_scores = collect(keys(element_scores))
    sort!(sorted_scores; by = x -> element_scores[x])

    for (a, b) = sorted_scores
        bsub, bpos = graph.mat_subposmap[b]
        if length(redundant_cols[(a, bsub)]) > 1 || length(redundant_rows[b]) > 1
            deleteat!(graph.clusters[a], findall(x->x==b,graph.clusters[a]))
            delete!(redundant_cols[(a, bsub)], (a, b))
            delete!(redundant_rows[b], (a, b))
        end
    end

    problem_cols = [(a, b) for (a, b) in keys(redundant_cols) if length(redundant_cols[(a, b)]) > 1]
    problem_rows = [a for a in keys(redundant_rows) if length(redundant_rows[a]) > 1]
    println("Finished. Now $(length(problem_rows)) row violations and $(length(problem_cols)) column violations")

    graph.clusters = [cluster for cluster in graph.clusters if length(cluster) > 1]
    println("Purged cluster violations. Found $(length(graph.clusters)) clean clusters..")
end

# https://discourse.julialang.org/t/minimum-by-key-function/56056/3
minimumby(f, iter) = reduce(iter) do x, y
    f(x) < f(y) ? x : y
end

function min_clusters_search(g :: AlnGraph)
    graph = g
    subset_clusters = TDict{Int, Vector{Tuple{Int, Int}}}()
    cluster_positions = TDict{Int, TDict{Int, Int}}()
    queue_idxs = TDict{Int, Int}()
    cluster_breaks = TDict{Tuple{Int, Int}, Vector{Int}}()
    max_frontier = TDict{Int, Int}()
    visited_states = Set{NTuple{25, Int64}}()
    maximal_cut = TDict{Int, Int}()
    state_counter = 0
    aggression = 1.0
    greedy = false
    totalpairs = 0
    last_frontier_state = State(
        ((0,0,0,0), 0, length(graph.clusters), totalpairs, state_counter, queue_idxs, cluster_breaks, maximal_cut, [], true))

    for (a, cluster) = enumerate(graph.clusters)
        for b = cluster
            bsub, bpos = graph.mat_subposmap[b]
            subset_clusters[bsub] = vcat(get(subset_clusters, bsub, Vector{Tuple{Int, Int}}()), [(a, bpos)])
            cluster_positions[a] = TDict()
            totalpairs += length(cluster) * (length(cluster) - 1) ÷ 2
        end
    end

    for asub = keys(subset_clusters)
        sort!(subset_clusters[asub]; by = x -> x[2])
        queue_idxs[asub] = 1

        max_frontier[asub] = -1
        maximal_cut[asub] = -1

        for i = 1:length(subset_clusters[asub])
            a = subset_clusters[asub][i][1]
            cluster_positions[a][asub] = i
        end
    end

    heap = BinaryMinHeap{State}()
    start_state = State(
        ((0,0,0,0), 0, length(graph.clusters), totalpairs, state_counter, queue_idxs, cluster_breaks, maximal_cut, [], true))
    start_state = develop_state(
        start_state, g, aggression, greedy, 0, subset_clusters, cluster_positions)
    push!(heap, start_state)

    while length(heap) > 0
        heapcleared = false
        if length(heap) > 5000
            if aggression == 1.0
                aggression = 1.2
                println("Increasing aggression to $(aggression)")
            elseif aggression < 8.0
                aggression = floor(aggression) * 2.0
                println("Increasing aggression to $(aggression)")
            else
                println("now full greedy")
                greedy = true
                aggression = 1.0
            end

            empty!(heap)
            empty!(visited_states)
            last_frontier_state = develop_state(last_frontier_state, graph, aggression, greedy, 0, subset_clusters, cluster_positions)
            push!(heap, last_frontier_state)
            heapcleared = true
        end

        state = pop!(heap)
        heuristic, numordered, numleft, pairsleft, counter, queue_idxs, cluster_breaks, maximal_cut, newcluster_breaks, safefrontier = state

        if isempty(newcluster_breaks)
            break
        else
            statekey = ([queue_idxs[sub] for sub = keys(subset_clusters)]...,)
            if statekey in visited_states
                continue
            else
                push!(visited_states, statekey)
            end

            new_sub_frontier = true
            for asub = keys(queue_idxs)
                if queue_idxs[asub] <= max_frontier[asub]
                    new_sub_frontier = false
                    break
                end
            end

            if new_sub_frontier
                max_frontier = queue_idxs
                println("Reached new search frontier")
                println(length(max_frontier))
                last_frontier_state = state
                greedy = false
            end

            if safefrontier &&! heapcleared
                println("Safe frontier reached.. Dumping $(length(heap)) from heap and resetting aggression..")
                last_frontier_state = state
                empty!(heap)
                empty!(visited_states)
                aggression = 1.0
                greedy = false
            end

            nextstates = []
            
            for (a, goodside, badside, crossedclusters) = newcluster_breaks
                g, b = length(goodside), length(badside)
                pairsdiff = g * (g - 1) / 2 + b * (b - 1) / 2 - (g+b)*(g+b-1)/2

                state_counter = state_counter + 1
                queue_idxs_copy = copy(queue_idxs)
                cluster_breaks_copy = copy(cluster_breaks)
                maximal_cut_copy = copy(maximal_cut)
                for b = goodside
                    bsub, bpos = graph.mat_subposmap[b]
                    cluster_breaks_copy[(a, bsub)] = goodside
                    maximal_cut_copy[bsub] = max(maximal_cut_copy[bsub], cluster_positions[a][bsub])
                end

                for b = badside
                    bsub, bpos = graph.mat_subposmap[b]
                    cluster_breaks_copy[(a, bsub)] = badside
                    maximal_cut_copy[bsub] = max(maximal_cut_copy[bsub], cluster_positions[a][bsub])
                end

                nextstate = State(
                    ((0,0,0,0), numordered, numleft+1, pairsleft + pairsdiff, state_counter, queue_idxs_copy, cluster_breaks_copy, maximal_cut_copy, [], false)
                )
                nextstate = develop_state(nextstate, graph, aggression, greedy, length(crossedclusters), subset_clusters, cluster_positions)

                push!(nextstates, nextstate)
            end

            if greedy
                nextstate = minimumby(x -> x[1], nextstates)
                push!(heap, nextstate)
            else
                for nextstate = nextstates
                    push!(heap, nextstate)
                end
            end
        end
    end

    queue_idxs = TDict()
    for asub = keys(subset_clusters)
        queue_idxs[asub] = 1
    end

    orderedclusters = []
    foundgood = true
    # @show sum(queue_idxs)
    # @show queue_idxs
    while foundgood
        foundgood = false
        t = length(orderedclusters)
        for asub = keys(queue_idxs)
            good = true
            idx = queue_idxs[asub]
            if idx > length(subset_clusters[asub])
                continue
            end
            a, pos = subset_clusters[asub][idx]

            if (a, asub) in keys(cluster_breaks)
                cluster = cluster_breaks[(a, asub)]
            else
                cluster = graph.clusters[a]
            end

            for b in cluster
                bsub, bpos = graph.mat_subposmap[b]
                if cluster_positions[a][bsub] != queue_idxs[bsub]
                    good = false
                    break
                end
            end
            
            if good
                push!(orderedclusters, cluster)
                for b in cluster
                    bsub, bpos = graph.mat_subposmap[b]
                    queue_idxs[bsub] = cluster_positions[a][bsub] + 1
                end
                # @show queue_idxs
                foundgood = true
                break
            end
        end
    end
    graph.clusters = orderedclusters
    return orderedclusters
end

function develop_state(
    state :: State, 
    graph :: AlnGraph, aggression :: Float64,
    greedy :: Bool, crossed :: Int, subset_clusters :: TDict{Int, Vector{Tuple{Int, Int}}}, cluster_positions :: TDict{Int, TDict{Int, Int}})
    heuristic, numordered, numleft, pairsleft, counter, queue_idxs, cluster_breaks, maximal_cut, new_cluster_breaks, safefrontier = state

    foundgood = true
    while foundgood
        foundgood = false
        new_cluster_breaks = Tuple{Int, Vector{Int}, Vector{Int}, Set{Int}}[]
        visited = Set{Tuple{Int, Int}}()
        safefrontier = true

        for asub = keys(queue_idxs)
            idx = queue_idxs[asub]
            if idx <= maximal_cut[asub]
                safefrontier = false
            end

            if idx > length(subset_clusters[asub])
                continue
            end
            a, pos = subset_clusters[asub][idx]
            if (a, asub) in visited
                continue
            end
            if (a, asub) in keys(cluster_breaks)
                cluster = cluster_breaks[(a, asub)]
            else
                cluster = graph.clusters[a]
            end

            goodside, badside, crossedclusters = Int[], Int[], Set{Int}()
            diffs = []
            for b in cluster
                bsub, bpos = graph.mat_subposmap[b]
                push!(visited, (a, bsub))
                bidx = cluster_positions[a][bsub]
                diff = bidx - queue_idxs[bsub]
                push!(diffs, diff)
                if diff == 0
                    push!(goodside, b)
                else
                    push!(badside, b)
                end
            end

            # @show diffs

            if length(badside) == 0
                for b = cluster
                    bsub, bpos = graph.mat_subposmap[b]
                    queue_idxs[bsub] = cluster_positions[a][bsub] + 1
                end
                numordered += 1
                numleft -= 1
                foundgood = true
                break
            else
                push!(new_cluster_breaks, (a, goodside, badside, crossedclusters))
            end
        end
    end

    if greedy
        for (a, goodside, badside, crossedclusters) = new_cluster_breaks
            goodsub = Set{Int}()
            for b = goodside
                bsub, bpos = graph.mat_subposmap[b]
                push!(goodsub, bsub)
            end

            for b = badside
                bsub, bpos = graph.mat_subposmap[b]
                for i = queue_idxs[bsub]:(cluster_positions[a][bsub])
                    c, posc = subset_clusters[bsub][i]
                    if (c, bsub) in keys(cluster_breaks)
                        othercluster = cluster_breaks[(c, bsub)]
                    else
                        othercluster = graph.clusters[c]
                    end

                    for csite = othercluster
                        csub, cpos = graph.mat_subposmap[csite]
                        if csub in goodsub && cluster_positions[c][csub] > cluster_positions[a][csub]
                            push!(crossedclusters, c)
                            break
                        end
                    end
                end
            end
        end
    end

    if safefrontier || length(new_cluster_breaks) == 0
        heuristic = (numleft + numordered, -numordered, -crossed, -pairsleft)
    else
        heuristic = (aggression * numleft + numordered, -numordered, -crossed, -pairsleft)
    end
    state = State(
        (heuristic, numordered, numleft, pairsleft,
        counter, queue_idxs, cluster_breaks, maximal_cut, new_cluster_breaks, safefrontier))
    return state
end

function dump_clusters_to_file(g :: AlnGraph, filename)
    open(filename, "w+") do f
        for cluster in g.clusters
            println(f, join([c - 1 for c in cluster], " "))
        end
    end
end

function convert_clusters_zerobased(g :: AlnGraph)
    for cluster in g.clusters
        cluster .-= 1
    end
    return g.clusters
end

function find_trace(c)
    g = AlnGraph(AlnContext(c.workingDir, c.subalignmentPaths))
    purge_duplicate_clusters(g)
    purge_cluster_violations(g)
    min_clusters_search(g)
    return g
end

precompile(AlnGraph, (AlnContext, ))
precompile(min_clusters_search, (AlnGraph,))
precompile(develop_state, 
  (State, AlnGraph, Float64, Int, TDict{Int, Vector{Tuple{Int, Int}}}, TDict{Int, TDict{Int, Int}}))

export min_clusters_search, develop_state, dump_clusters_to_file, AlnContext, AlnGraph
export purge_duplicate_clusters, purge_cluster_violations, convert_clusters_zerobased, find_trace
# using Traceur

# workingdir = "allofthem"
# using Glob
# b = glob(joinpath(workingdir, "subalignments","*"))
# sort!(b; by = x -> parse(Int, split(splitext(basename(x))[1], "_")[end]))
# context = AlnContext(workingdir, b)
# g = AlnGraph(context)
# purge_duplicate_clusters(g)
# purge_cluster_violations(g)
# min_clusters_search(g)
# dump_clusters_to_file(g, "clusters.scratch.txt")
end