# Steady-state timings (min of 3 runs; legacy: 1 run) for representative cases.
include(joinpath(@__DIR__, "..", "SpiderBench.jl")); using .SpiderBench
using Printf
run_new(1//7); run_new(1//7; feet=BigFloat); run_legacy(1//7)
out = open(joinpath(@__DIR__, "..", "results", "timing.csv"), "w")
println(out, "family,period,impl,status,iterations,min_time,bytes,points")
for fam in (:cusp, :tip, :approx16, :approx14, :doubling), p in (8, 16, 32, 64)
    fam === :doubling && !ispow2(p) && continue
    θ = family(fam, p)
    for (impl, f, reps) in (("legacy", () -> run_legacy(θ), 1),
                            ("new-f64", () -> run_new(θ; impl="new-f64"), 3),
                            ("new-bigfeet", () -> run_new(θ; impl="new-bigfeet", feet=BigFloat), 3))
        impl == "legacy" && p > 32 && continue
        rs = [f() for _ in 1:reps]
        r = rs[argmin([x.time for x in rs])]
        println(out, join([fam, p, impl, r.status, r.iterations, r.time, r.bytes, r.points], ","))
        @printf("%-9s p=%2d %-12s %-11s it=%4d  %.4f s  %8.1f MB  pts=%d\n", fam, p, impl, r.status, r.iterations, r.time, r.bytes/2^20, r.points)
        flush(stdout)
    end
end
close(out)
