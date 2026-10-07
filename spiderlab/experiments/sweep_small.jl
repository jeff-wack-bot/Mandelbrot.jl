# Experiment 1: exhaustive sweep over small periods / denominators.
# Usage: julia --project=spiderlab spiderlab/experiments/sweep_small.jl [pmax] [bmax]
include(joinpath(@__DIR__, "..", "SpiderBench.jl")); using .SpiderBench
pmax = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
bmax = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 64
angles = vcat([periodic_angles(p) for p in 2:pmax]..., preperiodic_angles(bmax))
println("$(length(angles)) angles")
# warm up
for f in (θ -> run_legacy(θ), θ -> run_legacy(θ; fixindex=true), θ -> run_new(θ)); f(1//7); f(1//6); end
rows = Row[]
t0 = time()
for (k, θ) in enumerate(angles)
    ref = reference(θ)
    ref.ok || println("reference not ok for $θ")
    push!(rows, run_legacy(θ))
    push!(rows, run_legacy(θ; fixindex=true))
    push!(rows, run_new(θ))
    push!(rows, run_new(θ; impl="new-nopolish", polish=false))
    k % 100 == 0 && (println("$k / $(length(angles))  $(round(time()-t0))s"); flush(stdout))
end
write_csv(joinpath(@__DIR__, "..", "results", "sweep_small.csv"), rows)
println("done in $(round(time()-t0))s")
