#!/bin/bash
#=
exec julia --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#
using PyCall
using Base.Filesystem
using Logging
pushfirst!(PyVector(pyimport("sys")["path"]), @__DIR__)
function main()
    external_tools = pyimport("tools.external_tools")
    # argument parsing, extremely ugly
    output_path = nothing
    output_ix = -1
    dir_ix = -1
    dir_path = nothing
    magus_args = copy(ARGS)
    its_ix = -1
    its_times = -1
    master_mode = "MAGUS_MASTER" ∈ keys(ENV)
    if master_mode
        info!("Running in master mode. Will do the FastTree things myself.")
    end
    for (j, i) = enumerate(magus_args)
        if i == "-o"
            output_ix = j
            output_path = magus_args[j+1]
        end
        if i == "-d"
            dir_ix = j
            dir_path = magus_args[j+1]
        end
        if i == "--its"
            its_ix = j
            its_times = parse(Int, magus_args[j+1])
        end
    end
    deleteat!(magus_args, sort([output_ix, output_ix + 1,dir_ix, dir_ix + 1, its_ix, its_ix + 1]))
    @info "Running MAGUS with total iterations: $(its_times)"
    @assert its_times > 1
    # we now have the output path, and the rest of the arguments.
    for i = 1:its_times
        @info "Starting iteration $i"
        output_filename = "$(output_path).it$(i)"
        output_treename = "$(output_path).it$(i).tre"
        env_dir = "$(dir_path)_it$(i)"
        initial_tree_arg = i == 1 ? [] : ["-t", "$(output_path).it$(i-1).tre"]
        if i <= 4
            run(`fmagus $magus_args -o $(output_filename) -d $(env_dir) $initial_tree_arg --graphclustermethod upgma --graphtracemethod None`)
        else
            run(`fmagus $magus_args -o $(output_filename) -d $(env_dir) $initial_tree_arg`)
        end
        if i < its_times # if we are not at the last iteration
            # we estimate the tree
            if master_mode
                @info "Compressing alignment: $(output_filename) -> $(output_filename).compressed"
                run(`ogcat mask -p 0.999 $output_filename -o $(output_filename).compressed`)
                task = external_tools.runFastTree(output_filename * ".compressed", env_dir, output_treename, "fast")
                taskArgs = task.taskArgs
                taskArgs["workingDir"] = "."
                task.taskArgs = taskArgs
                task.run()
            else
                @info "Waiting for FastTree to finish on master"
                while !isfile(output_treename)
                    @info "Waiting..."
                    sleep(5)
                end
                @info "Joining later..."
                sleep(8)
            end
        else
            if master_mode
                # last iteration, move the output file to the output path
                cp(output_filename, output_path)
            end
        end
        @info "Iteration $i finished."
    end
end

main()
