using Base.Filesystem
using Glob
using Distributed

# pattern = ARGS[1]
# @assert !isempty(pattern)

# WORKERS = 15
# if Sys.isapple()
#     WORKERS = 1
# end

# addprocs(WORKERS)

@everywhere TRUE_ALN_FILENAME = "true.indel.fa"
# @everywhere TRUE_ALN_FILENAME = "truealign.hf.txt"

@everywhere function parse_result(io :: String, keys)
    res = Dict()
    for l = eachline(IOBuffer(io))
        for k = keys
            if contains(l, k)
                # @show "contains!"
                src = split(l)[2]
                r = tryparse(Float64, src)
                res[k] = isnothing(r) ? src : r
            end
        end
    end
    return res
end

@everywhere function obtain_rates(ref, estimated)
    fastsp_script = joinpath(@__DIR__, "restricted_fastsp.py")
    res = read(`python3 $fastsp_script -r $ref -e $estimated`, String)
    r = parse_result(res, ["SPFP", "SPFN"])
    # @show res
    return (spfp = r["SPFP"], spfn = r["SPFN"])
end

tasks = []

@enum AlnType subset=1 support=2

struct Task
    envpath
    ref
    estimated
    type
end

@everywhere function execute_task(task :: Task)
    rates = obtain_rates(task.ref, task.estimated)
    return [task.estimated, task.type, (rates.spfp + rates.spfn)/2, rates.spfp, rates.spfn]
end



tasks = []
for arg in ARGS
    for p in glob(arg)
        graph_path = joinpath(p, "graph")
        subset_path = joinpath(p, "subalignments")
        subset_alignment_pattern = joinpath(subset_path, "subalignment_subset_*.txt")
        graph_alignment_pattern = joinpath(graph_path, "backbone_*_mafft.txt")
        tokens = splitpath(p)
        dataset = tokens[2]
        rep = tokens[3]
        true_aln_path = joinpath(dataset, rep, TRUE_ALN_FILENAME)
        for sa = glob(subset_alignment_pattern)
            push!(tasks, Task(p, true_aln_path, sa, 1))
        end
        for ga = glob(graph_alignment_pattern)
            push!(tasks, Task(p, true_aln_path, ga, 2))
        end
    end
end

@show tasks
# @show isempty(tasks)

results = map(execute_task, tasks)
println("inputpath,type,ref,error,spfp,spfn")
for r in results
    println(join(r, ","))
end