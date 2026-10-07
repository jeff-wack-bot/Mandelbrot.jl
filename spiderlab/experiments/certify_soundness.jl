# Experiment: soundness of `certify`. For every angle of period p and every nucleus of
# exact period p, ask certify(S_k, nucleus) for spiders S_k at several iterations k.
# A sound certificate accepts only the correct nucleus.
include(joinpath(@__DIR__, "..", "SpiderBench.jl")); using .SpiderBench
const M = SpiderBench.M
const VARIANTS = [("room+feet", (drag=:room, legs=false)), ("room+legs", (drag=:room, legs=true)),
                  ("isotopy+feet", (drag=:isotopy, legs=false)), ("isotopy+legs", (drag=:isotopy, legs=true))]
function main(pmax)
    out = open(joinpath(@__DIR__, "..", "results", "certify_soundness.csv"), "w")
    println(out, "variant,period,iteration,steps,candidates,true_certified,true_rejected,true_inconclusive,false_certified,false_rejected,false_inconclusive")
    for p in 2:pmax
        angles = periodic_angles(p)
        nuclei = ComplexF64[]
        for θ in angles
            c = ComplexF64(reference(θ).c)
            any(abs(c - d) < 1e-9 for d in nuclei) || push!(nuclei, c)
        end
        nuclei = vcat(nuclei, conj.(filter(c -> abs(imag(c)) > 1e-12, nuclei)))
        for (vname, vkw) in VARIANTS, steps in (1,), k in (1, 2, 3, 5, 10, 20)
            cnt = zeros(Int, 6)
            for θ in angles
                ctrue = ComplexF64(reference(θ).c)
                S = M.standardspider(θ)
                for _ in 1:k; S = M.spider_map(S); end
                for c in nuclei
                    v = M.certify(S, c; steps=steps, vkw...)
                    good = abs(c - ctrue) < 1e-9
                    idx = (good ? 0 : 3) + (v === :certified ? 1 : v === :rejected ? 2 : 3)
                    cnt[idx] += 1
                end
            end
            println(out, join([vname, p, k, steps, length(nuclei), cnt...], ","))
            println("$vname p=$p k=$k candidates=$(length(nuclei)) true(c/r/i)=$(cnt[1:3]) false(c/r/i)=$(cnt[4:6])"); flush(stdout)
        end
    end
    close(out)
end
main(length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 8)
