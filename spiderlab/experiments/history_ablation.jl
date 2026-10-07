# Experiment: (1) error and point-count histories, legacy vs upgraded spider map
#             (2) ablation of refinement / pruning on the exhaustive small set
include(joinpath(@__DIR__, "..", "SpiderBench.jl")); using .SpiderBench
const M = SpiderBench.M
const L = SpiderBench.LegacySpider
function histories()
    out = open(joinpath(@__DIR__, "..", "results", "history.csv"), "w")
    println(out, "angle,impl,iteration,err,points")
    L.FIX_GETINDEX[] = true
    for θ in (1//7, 9//56, 145//1023, 1//1023, family(:tip, 12), 1//4)
        ref = reference(θ)
        S = L.standardspider(Rational{Int}(θ))
        for k in 1:400
            S = L.spider_map(S)
            println(out, "$θ,legacy+fix,$k,$(Float64(abs(S.legs[2][end]/2 - ref.c))),$(sum(length, S.legs))")
        end
        for (name, kw) in (("new", (;)), ("new-noprune", (prune=false,)), ("new-norefine", (refine=false,)))
            S = M.standardspider(θ)
            for k in 1:400
                S = try M.spider_map(S; maxpoints=2_000_000, kw...) catch; break end
                e = Float64(abs(S.feet[2]/2 - ref.c))
                println(out, "$θ,$name,$k,$e,$(sum(length, S.legs))")
                e < 1e-15 && break
            end
        end
        println("history $θ done"); flush(stdout)
    end
    close(out)
end
function ablation()
    angles = vcat([periodic_angles(p) for p in 2:9]..., preperiodic_angles(48))
    rows = Row[]
    for θ in angles
        for (name, kw) in (("new", (;)), ("new-noprune", (prune=false,)), ("new-norefine", (refine=false,)),
                           ("new-neither", (prune=false, refine=false)))
            push!(rows, run_new(θ; impl=name, maxpoints=2_000_000, kw...))
            push!(rows, run_new(θ; impl=name * "-nopolish", polish=false, maxpoints=2_000_000, kw...))
        end
    end
    write_csv(joinpath(@__DIR__, "..", "results", "ablation.csv"), rows)
    println("ablation done")
end
histories()
ablation()
