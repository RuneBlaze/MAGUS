using Profile
push!(LOAD_PATH, "MagusNight/")
using MagusNight
using Glob
using Profile
using BenchmarkTools
using Traceur

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

function main_task()
    workingdir = "../../Downloads/for_chengze/magus10krun_norecurse/"
    b = glob(joinpath(workingdir, "subalignments","*"))
    sort!(b; by = x -> parse(Int, split(splitext(basename(x))[1], "_")[end]))
    context = AlnContext(workingdir, b)
    g = AlnGraph(context)
    purge_duplicate_clusters(g)
    purge_cluster_violations(g)
    min_clusters_search(g)
    dump_clusters_to_file(g, "clusters.julia.txt")
end

run_benchmark()
# Profile.print()