# a hack to run GCM
using Logging
using Base.Filesystem
using Glob

function run_mafft(inpath, outpath)
    run(
        pipeline(`mafft --localpair --maxiterate 1000 --ep 0.123 --quiet --thread 16 --anysymbol $inpath`, outpath))
end

unaligned_filename = ARGS[1]
aligned_filename = ARGS[2]
env_path = ARGS[1] * "_env"
@info "Running GCM on $(unaligned_filename) in $(env_path) => $(aligned_filename)"
subtree_path = unaligned_filename * ".tre"
@assert isfile(subtree_path) "Subtree file $(subtree_path) does not exist"
run(`gcm137 slice -g 10x100 -s 100 -i $unaligned_filename -t $subtree_path -o $env_path`)
for unaligned = glob(joinpath(env_path, "*", "*.unaln.fa"))
    d, bn = splitdir(unaligned)
    nbn = replace(bn, ".unaln.fa" => "aln.fa")
    @info "Running MAFFT on $(unaligned) -> $(joinpath(d, nbn))"
    opath = joinpath(d, nbn)
    run_mafft(unaligned, opath)
end
@info "MAFFT-linsi runs finished"
aligned_constraints = glob(joinpath(env_path, "constraints", "*.aln.fa"))
aligned_glues = glob(joinpath(env_path, "glues", "*.aln.fa"))
run(`gcm137 merge -i $aligned_constraints -g $aligned_glues -o $aligned_filename`)