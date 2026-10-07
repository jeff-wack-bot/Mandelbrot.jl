# Profiling / micro-benchmarks of the upgraded spider map.
using Mandelbrot, Profile, BenchmarkTools, Printf
const M = Mandelbrot
include(joinpath(@__DIR__, "..", "SpiderBench.jl")); using .SpiderBench: family
function spider_after(θ, k; kw...)
    S = M.standardspider(θ; kw...)
    for _ in 1:k; S = M.spider_map(S); end
    S
end
for (name, θ, kw) in (("cusp p=30, F64", family(:cusp, 30), (;)),
                      ("tip p=30, BigFloat feet", family(:tip, 30), (feet=BigFloat,)),
                      ("tip p=30, BigFloat feet+legs", family(:tip, 30), (feet=BigFloat, legs=BigFloat)))
    S = spider_after(θ, 20; kw...)
    b = @benchmark M.spider_map($S) samples=50 evals=1
    @printf("%-32s points=%5d  spider_map median %.3f ms, %.1f KiB, %d allocs\n", name, sum(length, S.legs),
            1e3 * median(b.times) / 1e9 * 1e3 / 1e3, median(b).memory / 1024, median(b).allocs)
end
# where does the time go? (cusp p=30 without Newton: the long-run case)
θ = family(:cusp, 30)
M.spider_solve(θ; polish=false, maxiter=10)
Profile.clear()
@profile M.spider_solve(θ; polish=false, maxiter=1500)
Profile.print(format=:flat, sortedby=:count, mincount=40, noisefloor=2)
