#!/bin/bash
#=
exec julia --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#
numworkers = Sys.CPU_THREADS ÷ 4 + 1 # we overload the CPU by 1 more worker
cd(@__DIR__) do
    cmd = `huey_consumer.py tasks.remote_tasks.huey -n -w $numworkers -k process`
    run(cmd)
end