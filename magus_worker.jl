#!/bin/bash
#=
exec julia --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#
numworkers = Sys.CPU_THREADS ÷ 4 + 1 # we overload the CPU by 1 more worker
function run_with_timeout(command, timeout::Integer = 5)
    cmd = run(command; wait=false)
    for _ in 1:timeout
        if !process_running(cmd) return success(cmd) end
        sleep(1)
    end
    kill(cmd)
    return false
end
cd(@__DIR__) do
    cmd = `huey_consumer.py tasks.remote_tasks.huey -n -w $numworkers -k process`
    run_with_timeout(cmd, trunc(Int, 3600 * 3.9 + rand(Float64) * 30))
end