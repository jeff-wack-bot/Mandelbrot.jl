# Experiment: scaling with the period for structured angle families.
# Usage: julia --project=spiderlab spiderlab/experiments/scaling.jl
include(joinpath(@__DIR__, "..", "SpiderBench.jl")); using .SpiderBench
using Random
const FAMILIES = (:cusp, :tip, :approx16, :approx14, :random, :doubling)
const PERIODS = [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]
const LEGACY_MAX = 32      # legacy cost grows quickly; beyond this it is skipped
function main()
    rows = Row[]
    run_new(1//7); run_new(1//7; feet=BigFloat); run_legacy(1//7)   # warm up
    rng = MersenneTwister(1)
    for fam in FAMILIES, p in PERIODS
        fam === :doubling && !ispow2(p) && continue
        θs = fam === :random ? [family(:random, p; rng=rng) for _ in 1:5] : [family(fam, p)]
        for θ in θs
            ref = reference(θ)
            tag = string(fam)
            p <= LEGACY_MAX && (r = run_legacy(θ); r.note = tag * " " * r.note; push!(rows, r))
            r = run_new(θ; impl="new-f64"); r.note = tag * " " * r.note; push!(rows, r)
            r = run_new(θ; impl="new-bigfeet", feet=BigFloat); r.note = tag * " " * r.note; push!(rows, r)
            if p <= 32
                r = run_new(θ; impl="new-bigall", feet=BigFloat, legs=BigFloat); r.note = tag * " " * r.note; push!(rows, r)
            end
            println(fam, " ", p, " ", θ, " size=", ref.size, " ",
                    join(["$(x.impl):$(x.status)/$(x.iterations)/$(round(x.time,digits=3))s" for x in rows[end-2:end]], "  "))
            flush(stdout)
        end
    end
    write_csv(joinpath(@__DIR__, "..", "results", "scaling.csv"), rows)
end
main()
