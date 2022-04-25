using JobSchedulers, Dates, Pipelines
using Logging
using DataFrames

set_scheduler_max_mem()
scheduler_start()
while scheduler_status(verbose=false) == :not_running
    @info "waiting for scheduler to properly start"
    sleep(1)
end
scheduler_status()
if Sys.islinux()
    C = Threads.nthreads() - 1
    tentative = C ÷ 3
    if mod(C, 4) < mod(C, 3)
        tentative = C ÷ 4
    end
    MAFFT_NUMTHREADS = min(tentative, 10)
else
    MAFFT_NUMTHREADS = 1 # for debugging, assert that each mafft process will only take 1 thread
end

@info "Assigned MAFFT_NUMTHREADS = $MAFFT_NUMTHREADS"

@info "Started internal job queue for MAFFT..."
function request_mafft(frompath, topath; name = "something", priority = 2)
    p = CmdProgram(
        inputs = "INPATH",
        outputs = "OUTPATH",
        cmd = pipeline(`mafft --localpair --maxiterate 1000 --ep 0.123 --quiet --thread $MAFFT_NUMTHREADS --anysymbol INPATH`, "OUTPATH"),
    )
    inputs = "INPATH" => frompath
    outputs = "OUTPATH" => topath
    job = Job(p, inputs, outputs; name = name, ncpu = MAFFT_NUMTHREADS, priority = priority, touch_run_id_file = false)
    submit!(job)
end

function request_constraint(frompath, topath)
    # scheduler_start()
    request_mafft(frompath, topath; name = "constraint", priority = 1)
end

function request_backbone(frompath, topath)
    # scheduler_start()
    request_mafft(frompath, topath; name = "backbone", priority = 2)
end

function check_constraints_finished()
    # @show queue()
    
    return isempty(filter(:name => n -> n == "constraint", queue()))
end

function check_backbones_finished()
    # @show queue()
    return isempty(filter(:name => n -> n == "constraint", queue()))
end

function init_scheduler()
    scheduler_start()
end

function everything_finished()
    # update_queue!()
    stats = nrow(queue())
    @info stats
    return stats == 0
end

function wait_for_everything()
    while !everything_finished()
        # scheduler_wait()
        JobSchedulers.update_queue!()
        sleep(3)
    end
end

# if abspath(PROGRAM_FILE) == @__FILE__
#     request_constraint("example/small.txt", "example/small_constraint.txt")
#     while !check_constraints_finished()
#         sleep(1)
#     end
# end