# Experiment: soundness of `certify` (default settings) against harder candidate sets:
#  (a) periodic angles, candidates = nuclei of every period dividing p
#  (b) preperiodic angles (denominator <= bmax), candidates = every Misiurewicz point with
#      preperiod <= l and period dividing k, plus every nucleus with period dividing k.
include(joinpath(@__DIR__, "..", "SpiderBench.jl")); using .SpiderBench
const M = SpiderBench.M
const GT = SpiderBench.GroundTruth
addpt!(v, c) = any(abs(c - d) < 1e-9 for d in v) || push!(v, c)
function nuclei_dividing(p)
    v = ComplexF64[]
    for d in 1:p
        p % d == 0 || continue
        d == 1 && (addpt!(v, 0.0im); continue)
        for θ in periodic_angles(d)
            c = ComplexF64(reference(θ).c); addpt!(v, c); addpt!(v, conj(c))
        end
    end
    return v
end
function main(pmax, bmax)
    out = open(joinpath(@__DIR__, "..", "results", "certify_soundness2.csv"), "w")
    println(out, "kind,type,iteration,angles,candidates,true_certified,true_rejected,true_inconclusive,false_certified,false_rejected,false_inconclusive")
    ks = (1, 2, 3, 5, 10, 20)
    for p in 4:pmax
        cands = nuclei_dividing(p)
        for k in ks
            cnt = zeros(Int, 6)
            for θ in periodic_angles(p)
                ctrue = ComplexF64(reference(θ).c)
                S = M.standardspider(θ); for _ in 1:k; S = M.spider_map(S); end
                for c in cands
                    v = M.certify(S, c); good = abs(c - ctrue) < 1e-9
                    cnt[(good ? 0 : 3) + (v === :certified ? 1 : v === :rejected ? 2 : 3)] += 1
                end
            end
            println(out, join(["periodic", p, k, length(periodic_angles(p)), length(cands), cnt...], ","))
            println("periodic p=$p k=$k cands=$(length(cands)) true=$(cnt[1:3]) false=$(cnt[4:6])"); flush(stdout)
        end
    end
    pre = preperiodic_angles(bmax)
    types = Dict{Tuple{Int,Int},Vector{Rational{BigInt}}}()
    for θ in pre; push!(get!(types, GT.angle_type(θ), Rational{BigInt}[]), θ); end
    mis = Dict{Tuple{Int,Int},Vector{ComplexF64}}()
    for (t, θs) in types, θ in θs
        v = get!(mis, t, ComplexF64[]); c = ComplexF64(reference(θ).c); addpt!(v, c); addpt!(v, conj(c))
    end
    for t in sort(collect(keys(types)))
        l, kper = t
        cands = ComplexF64[]
        for (t2, v) in mis
            (t2[1] <= l && kper % t2[2] == 0) || continue
            for c in v; addpt!(cands, c); end
        end
        for c in nuclei_dividing(kper); addpt!(cands, c); end
        for k in ks
            cnt = zeros(Int, 6)
            for θ in types[t]
                ctrue = ComplexF64(reference(θ).c)
                S = M.standardspider(θ); for _ in 1:k; S = M.spider_map(S); end
                for c in cands
                    v = M.certify(S, c); good = abs(c - ctrue) < 1e-9
                    cnt[(good ? 0 : 3) + (v === :certified ? 1 : v === :rejected ? 2 : 3)] += 1
                end
            end
            println(out, join(["preperiodic", "$(l)_$(kper)", k, length(types[t]), length(cands), cnt...], ","))
            println("preperiodic $t k=$k angles=$(length(types[t])) cands=$(length(cands)) true=$(cnt[1:3]) false=$(cnt[4:6])"); flush(stdout)
        end
    end
    close(out)
end
main(9, 64)
