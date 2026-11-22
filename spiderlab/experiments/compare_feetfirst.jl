# Feet-first (current src/spidermap.jl) vs shoulder-first (UpgradeSpider = 38caca9):
# (1) exhaustive sweep, (2) scaling families, (3) steady-state timings.
# Usage: julia --project=spiderlab spiderlab/experiments/compare_feetfirst.jl [sweep|scaling|timing]
include(joinpath(@__DIR__, "..", "SpiderBench.jl")); using .SpiderBench
using Mandelbrot
using Random, Printf
const OUT = joinpath(@__DIR__, "..", "results")
const IMPLS = (("upgrade", UpgradeSpider), ("feetfirst", Mandelbrot))
what = isempty(ARGS) ? "sweep" : ARGS[1]
for (_, m) in IMPLS; run_new(1//7; mod=m); run_new(1//6; mod=m); run_new(1//7; mod=m, feet=BigFloat); end

if what == "sweep"
    angles = vcat([periodic_angles(p) for p in 2:10]..., preperiodic_angles(64))
    rows = Row[]; t0 = time()
    for (k, θ) in enumerate(angles)
        reference(θ)
        for (name, m) in IMPLS
            push!(rows, run_new(θ; impl=name, mod=m))
            push!(rows, run_new(θ; impl=name * "-nopolish", mod=m, polish=false))
        end
        push!(rows, run_new(θ; impl="feetfirst-norefine-nopolish", refine=false, polish=false))
        k % 200 == 0 && (println("$k / $(length(angles))  $(round(time()-t0))s"); flush(stdout))
    end
    write_csv(joinpath(OUT, "compare_sweep.csv"), rows)
elseif what == "scaling"
    rows = Row[]; rng = MersenneTwister(1)
    for fam in (:cusp, :tip, :approx16, :approx14, :random, :doubling), p in [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]
        fam === :doubling && !ispow2(p) && continue
        θs = fam === :random ? [family(:random, p; rng=rng) for _ in 1:5] : [family(fam, p)]
        for θ in θs
            reference(θ)
            for (name, m) in IMPLS, (tag, kw) in (("f64", (;)), ("bigfeet", (feet=BigFloat,)))
                r = run_new(θ; impl="$name-$tag", mod=m, kw...); r.note = "$fam " * r.note; push!(rows, r)
            end
            println(fam, " ", p, " ", join(["$(x.impl):$(x.status)/$(x.iterations)" for x in rows[end-3:end]], "  ")); flush(stdout)
        end
    end
    write_csv(joinpath(OUT, "compare_scaling.csv"), rows)
elseif what == "timing"
    open(joinpath(OUT, "compare_timing.csv"), "w") do out
        println(out, "family,period,impl,status,iterations,min_time,bytes,points")
        for fam in (:cusp, :tip, :approx16, :approx14, :doubling), p in (8, 16, 32, 64)
            fam === :doubling && !ispow2(p) && continue
            θ = family(fam, p)
            for (name, m) in IMPLS, (tag, kw) in (("f64", (;)), ("bigfeet", (feet=BigFloat,)))
                rs = [run_new(θ; impl="$name-$tag", mod=m, kw...) for _ in 1:3]
                r = rs[argmin([x.time for x in rs])]
                println(out, join([fam, p, r.impl, r.status, r.iterations, r.time, r.bytes, r.points], ","))
                @printf("%-9s p=%2d %-18s %-11s it=%4d  %.4f s  %7.2f MB  pts=%d\n", fam, p, r.impl, r.status, r.iterations, r.time, r.bytes/2^20, r.points)
                flush(stdout)
            end
        end
    end
end
