# Per-step correctness of the lifting. For each angle, iterate each implementation; at every
# step compare its output with a reference step computed from the *same* input spider:
# legs densified (each piece shorter than 1e-3 of its distance to 0), then one step of the
# shoulder-first map with refinement and pruning switched off. A leg is counted as wrong if
# its winding numbers against the reference leg (around the other feet) are not all equal.
# Usage: julia --project=spiderlab spiderlab/experiments/lift_correctness.jl [steps [pmin:pmax [k/m]]]
# With a period range only periodic angles of those periods are used (part k of m).
include(joinpath(@__DIR__, "..", "SpiderBench.jl")); using .SpiderBench
using Mandelbrot
const M = Mandelbrot; const U = UpgradeSpider
steps = isempty(ARGS) ? 30 : parse(Int, ARGS[1])
prange = length(ARGS) >= 2 ? (:)(parse.(Int, split(ARGS[2], ":"))...) : nothing
part = length(ARGS) >= 3 ? parse.(Int, split(ARGS[3], "/")) : [1, 1]

function dense(P)
    out = [P[1]]
    function rec(a, b)
        if abs(b - a) > 1e-3 * min(abs(a), abs(b)) + 1e-9
            m = (a + b) / 2; rec(a, m); rec(m, b)
        else
            push!(out, b)
        end
    end
    for k in 1:length(P)-1; rec(P[k], P[k+1]); end
    return out
end
reference_step(S) = U.spider_map(U.Spider{Float64,Float64}(S.angle, S.orbit, S.kneading_sequence,
                                 [dense(P) for P in S.legs], copy(S.feet)); refine=false, prune=false)
function wronglegs(S, R)
    bad = 0
    for i in eachindex(S.feet)
        w = M.leg_windings(S.legs[i], R.legs[i], S.feet, i)
        (w === nothing || !allequal([w[j] for j in eachindex(w) if j != i])) && (bad += 1)
    end
    return bad
end

angles = prange === nothing ? vcat([periodic_angles(p) for p in 2:8]..., preperiodic_angles(32)) :
         vcat([periodic_angles(p) for p in prange]...)[part[1]:part[2]:end]
fname = prange === nothing ? "lift_correctness" : "lift_correctness_p$(replace(ARGS[2], ":" => "-"))_$(part[1])of$(part[2])"
out = open(joinpath(@__DIR__, "..", "results", fname * ".csv"), "w")
println(out, "impl,angle,steps,bad_steps,bad_legs,feet_mismatch,first_bad")
for θ in angles, (name, mod) in (("upgrade", U), ("feetfirst", M))
    S = mod.standardspider(Rational{Int}(θ))
    badsteps = 0; badlegs = 0; feetbad = 0; first = 0
    for it in 1:steps
        R = reference_step(S)
        S = mod.spider_map(S)
        fm = maximum(abs.(S.feet .- R.feet)) > 1e-8 * max(1, maximum(abs, R.feet))
        b = fm ? 0 : wronglegs(S, R)
        if fm || b > 0
            badsteps += 1; badlegs += b; feetbad += fm; first == 0 && (first = it)
        end
    end
    println(out, join([name, θ, steps, badsteps, badlegs, feetbad, first], ","))
    flush(out)
end
close(out)
println("done")
