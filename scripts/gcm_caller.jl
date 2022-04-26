using Distributed
using Base.Filesystem
using Logging
basedir = ARGS[1]
mstdir = joinpath(basedir, "decomposition", "mst.tre")

PAIRWISE_GEN = joinpath(@__DIR__, "pairwise_merger.py")
run(`python3 $PAIRWISE_GEN -i $mstdir`)

@info "Generated info files..."

merge_commands_file = joinpath(basedir, "merged", "merger.txt")
@assert isfile(merge_commands_file)

merge_order_file = joinpath(basedir, "merged", "order.txt")
@assert isfile(merge_order_file)

@info "Found merge order file... at $(merge_order_file)"

merge_commands = [] # list of four-tuple of strings
merge_order = [] # list of strings

open(merge_commands_file) do f
    for l in eachline(f)
        push!(merge_commands, split(l))
    end
end

open(merge_order_file) do f
    for l in eachline(f)
        push!(merge_order, strip(l))
    end
end

addprocs(10)

@everywhere function merge_aln(lhs, rhs, supp, output)
    run(`magus -np 1 --subalignments $lhs $rhs -b $supp -d $(output)_env -o $output`)
end

@info "Started pairwise merging"

@sync @distributed for command = merge_commands
    lhs = command[1]
    rhs = command[2]
    supp = command[3]
    output = command[4]
    merge_aln(lhs, rhs, supp, output)
end

@info "Finished pairwise merging. Generated type 2 subalignments."