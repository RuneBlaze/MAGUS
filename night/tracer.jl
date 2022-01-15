using Profile
push!(LOAD_PATH, "MagusNight/")
using MagusNight
using Glob
using Profile
using BenchmarkTools
using Traceur
using ProfileView

function run_benchmark()
    workingdir = "../../Downloads/sandia_data/magus10krun_norecurse/"
    b = glob(joinpath(workingdir, "subalignments","*"))
    sort!(b; by = x -> parse(Int, split(splitext(basename(x))[1], "_")[end]))
    context = AlnContext(workingdir, b)
    b = @benchmark begin
        g = AlnGraph($context)
        purge_duplicate_clusters(g)
        purge_cluster_violations(g)
        min_clusters_search(g)
    end seconds=60
    io = IOBuffer()
    show(io, "text/plain", b)
    s = String(take!(io))
    println(s)
end

function debug_print_trace()
    workingdir = "../../Downloads/sandia_data/magus10krun_norecurse/"
    b = glob(joinpath(workingdir, "subalignments","*"))
    sort!(b; by = x -> parse(Int, split(splitext(basename(x))[1], "_")[end]))
    context = AlnContext(workingdir, b)
    g = AlnGraph(context)
    open("scratch/trace.txt", "w") do f
        debug_print_mapped_graph(f, g)
    end
end

function find_clusterings()
    workingdir = "../../Downloads/sandia_data/magus10krun_norecurse/"
    b = glob(joinpath(workingdir, "subalignments","*"))
    sort!(b; by = x -> parse(Int, split(splitext(basename(x))[1], "_")[end]))
    context = AlnContext(workingdir, b)
    labels, adj = read_graph(get_graph_path(context))
    g = AlnGraph(context)
    outfile = "scratch/clustering_results.test.txt"
    if isfile(outfile) && false
        results = open(outfile, "r") do f
            read(f, String)
        end
        results = eval(Meta.parse(results))
    else
        results = upgma_naive_clustering(labels, adj, g; config = ClusteringConfig(true, true))
        open(outfile, "w+") do f
            println(f, results)
        end
    end
    
    flatclusters = convert_to_flatclusters(results, g)
    @show check_flatclusters_validity(flatclusters)
end

function find_fast_clusterings()
    workingdir = "../../Downloads/sandia_data/magus10krun_norecurse/"
    b = glob(joinpath(workingdir, "subalignments","*"))
    sort!(b; by = x -> parse(Int, split(splitext(basename(x))[1], "_")[end]))
    context = AlnContext(workingdir, b)
    labels, adj = read_graph(get_graph_path(context))
    g = AlnGraph(context)
    results = @time fast_upgma(labels, adj, g; config = ClusteringConfig(true, true))
    # open("scratch/profile.txt", "w+") do f
    #     Profile.print(f)
    # end
    # @time fast_upgma(labels, adj, g; config = ClusteringConfig(true, true))
    # io = IOBuffer()
    # show(io, "text/plain", bench)
    # println(String(take!(io)))
    flatclusters = convert_to_flatclusters(results, g)
    @show check_flatclusters_validity(flatclusters)
end

function find_decompositions()
    workingdir = "../../Downloads/sandia_data/magus10krun_norecurse/"
    b = glob(joinpath(workingdir, "subalignments","*"))
    sort!(b; by = x -> parse(Int, split(splitext(basename(x))[1], "_")[end]))
    context = AlnContext(workingdir, b)
    labels, adj = read_graph(get_graph_path(context))
    g = AlnGraph(context)
    cc = [(length(e.rows), length(e.nodes)) for e = connected_components(labels, adj, g)]
    println(cc)
end

function main_task()
    workingdir = "../../Downloads/sandia_data/magus10krun_norecurse/"
    b = glob(joinpath(workingdir, "subalignments","*"))
    sort!(b; by = x -> parse(Int, split(splitext(basename(x))[1], "_")[end]))
    context = AlnContext(workingdir, b)
    g = AlnGraph(context)
    purge_duplicate_clusters(g)
    purge_cluster_violations(g)
    min_clusters_search(g)
    dump_clusters_to_file(g, "scratch/clusters.julia.txt")
end

# b = @profile begin
# debug_print_trace()
# println("hello")
# s = readline()
find_fast_clusterings()
# end

# io = IOBuffer()
# show(io, "text/plain", b)
# s = String(take!(io))
# println(s)

# run_benchmark()
# with open("scratch/profile.txt", "w+") do f
#     Profile.print(f)
# end