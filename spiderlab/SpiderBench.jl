# Harness for testing / benchmarking the spider algorithm.
#
#   include("SpiderBench.jl"); using .SpiderBench
#
# Provides angle families, a uniform runner for the legacy and the upgraded
# implementation, and result classification against GroundTruth.
module SpiderBench

using Mandelbrot
using Printf
using Random

const M = Mandelbrot
include(joinpath(@__DIR__, "LegacySpider.jl"))
include(joinpath(@__DIR__, "UpgradeSpider.jl"))
include(joinpath(@__DIR__, "GroundTruth.jl"))
using .GroundTruth

export LegacySpider, UpgradeSpider, GroundTruth, run_legacy, run_new, classify, reference,
       periodic_angles, preperiodic_angles, family, Row, write_csv, read_csv

#---------------------------------------------------------------------------
# angle families

"All angles in (0, 1/2] of exact period p (denominator 2^p - 1)."
function periodic_angles(p::Int)
    d = big(2)^p - 1
    out = Rational{BigInt}[]
    for a in 1:(d÷2)
        θ = a // d
        GroundTruth.angle_type(θ) == (0, p) && push!(out, θ)
    end
    return out
end

"All strictly preperiodic angles a/b in (0, 1/2] with b ≤ bmax."
function preperiodic_angles(bmax::Int)
    out = Set{Rational{BigInt}}()
    for b in 2:bmax, a in 1:(b÷2)
        θ = big(a) // b
        GroundTruth.angle_type(θ)[1] > 0 && push!(out, θ)
    end
    return sort!(collect(out))
end

bits2angle(bits::AbstractString) = parse(BigInt, bits; base=2) // (big(2)^length(bits) - 1)

"""
    family(name, p; rng)

Periodic angles of period p from named families:
- :cusp      1/(2^p-1)                 (near the cusp 1/4, slow convergence)
- :tip       (2^(p-1)-1)/(2^p-1)       (real axis, next to c = -2; components ~ 4^-p)
- :approx16  first p bits of 1/6 = .0010101..., made periodic (near c = i)
- :approx14  .0100...0 (period p, near the Misiurewicz point of 1/4)
- :random    uniformly random angle of exact period p
- :doubling  period 2^k angles of the period doubling cascade (p must be 2^k)
"""
function family(name::Symbol, p::Int; rng=Random.default_rng())
    if name === :cusp
        return big(1) // (big(2)^p - 1)
    elseif name === :tip
        return (big(2)^(p - 1) - 1) // (big(2)^p - 1)
    elseif name === :approx16
        bits = "00" * repeat("10", p)
        return bits2angle(bits[1:p-1] * "1")
    elseif name === :approx14
        return bits2angle("01" * repeat("0", p - 2))
    elseif name === :random
        while true
            θ = rand(rng, 1:(big(2)^p - 2)) // (big(2)^p - 1)
            GroundTruth.angle_type(θ) == (0, p) && return θ
        end
    elseif name === :doubling
        # tuning: replace 0 -> 01, 1 -> 10 starting from "01" (1/3)
        w = "01"
        while length(w) < p
            w = join(ch == '0' ? "01" : "10" for ch in w)
        end
        length(w) == p || error("p must be a power of 2")
        return bits2angle(w)
    end
    error("unknown family $name")
end

#---------------------------------------------------------------------------
# reference values (cached)

const REFCACHE = Dict{Rational{BigInt},Any}()
function reference(θ::Rational; prec::Int=256)
    θ = Rational{BigInt}(θ)
    get!(REFCACHE, θ) do
        l, p = GroundTruth.angle_type(θ)
        pr = max(prec, 8 * (l + p) + 64)
        r = GroundTruth.reference_parameter(θ; prec=pr)
        size = l == 0 ? Float64(setprecision(() -> GroundTruth.size_estimate(r.c, p), BigFloat, pr)) : NaN
        (c=r.c, ok=r.ok, preperiod=l, period=p, size=size)
    end
end

#---------------------------------------------------------------------------
# results

Base.@kwdef mutable struct Row
    impl::String = ""
    angle::String = ""
    preperiod::Int = 0
    period::Int = 0
    status::String = ""          # ok / inaccurate / wrong / nan / nonconverged / exception
    err::Float64 = NaN           # |c - c_ref|
    relerr::Float64 = NaN        # err / size of the component (periodic) or err / max(1,|c|)
    iterations::Int = 0
    time::Float64 = NaN          # seconds
    bytes::Float64 = NaN         # bytes allocated
    points::Int = 0              # points in the final spider
    maxpoints::Int = 0
    finish::String = ""
    c_re::Float64 = NaN
    c_im::Float64 = NaN
    note::String = ""
end

const FIELDS = fieldnames(Row)

function write_csv(path, rows::Vector{Row})
    open(path, "w") do io
        println(io, join(FIELDS, ","))
        for r in rows
            println(io, join((replace(string(getfield(r, f)), "," => ";") for f in FIELDS), ","))
        end
    end
end

function read_csv(path)
    lines = readlines(path)
    hdr = Symbol.(split(lines[1], ","))
    rows = Row[]
    for ln in lines[2:end]
        vals = split(ln, ",")
        r = Row()
        for (h, v) in zip(hdr, vals)
            T = fieldtype(Row, h)
            setfield!(r, h, T === String ? String(v) : parse(T, v))
        end
        push!(rows, r)
    end
    return rows
end

"""
classify(c, ref): ok if |c - c_ref| <= tol_ok·max(1,|c|); wrong if the error is larger
than the component size (periodic) or 1e-3 (preperiodic), i.e. a different
component/point; inaccurate otherwise.
"""
function classify!(row::Row, c, ref; tol_ok=1e-10)
    if c === nothing || !isfinite(c)
        row.status = "nan"
        return row
    end
    e = Float64(abs(Complex{BigFloat}(c) - ref.c))
    row.err = e
    scale = max(1.0, Float64(abs(ref.c)))
    wrongtol = isnan(ref.size) ? 1e-3 : min(1e-3, ref.size / 4)
    row.relerr = isnan(ref.size) ? e / scale : e / ref.size
    # ok: accurate to tol_ok AND inside the right component; unresolved: as accurate as
    # Float64 allows but the component is smaller than that; wrong: another component
    row.status = e > wrongtol ? (e <= 1e3 * eps(Float64) * scale ? "unresolved" : "wrong") :
                 e <= tol_ok * scale ? "ok" : "inaccurate"
    row.c_re = Float64(real(c)); row.c_im = Float64(imag(c))
    return row
end

#---------------------------------------------------------------------------
# runners

"""
Legacy implementation exactly as `parameter(theta, maxiter)` at 619df42 (absolute
stopping tolerance 1e-15), instrumented. `fixindex=true` uses the corrected Sequence
indexing for the kneading symbols.
"""
function run_legacy(θ::Rational; maxiter::Int=250, fixindex::Bool=false, tol_ok=1e-10,
                    maxpoints::Int=20_000_000)
    ref = reference(θ)
    row = Row(impl=fixindex ? "legacy+fix" : "legacy", angle=string(θ), preperiod=ref.preperiod,
              period=ref.period)
    LegacySpider.FIX_GETINDEX[] = fixindex
    if ref.period > 62 || denominator(θ) > big(2)^62
        row.status = "exception"; row.note = "Rational{Int} overflow"; return row
    end
    θi = Rational{Int}(θ)
    local c, it, pts, maxpts
    GC.gc()
    stats = @timed begin
        S = LegacySpider.standardspider(θi)
        c_last = S.legs[2][end] / 2
        c = c_last
        it = 0
        maxpts = sum(length, S.legs)
        converged = false
        try
            for ii in 1:maxiter
                S = LegacySpider.spider_map(S)
                it = ii
                c = S.legs[2][end] / 2
                np = sum(length, S.legs)
                maxpts = max(maxpts, np)
                if abs(c - c_last) < 1e-15
                    converged = true
                    break
                end
                np > maxpoints && error("too many points")
                c_last = c
            end
        catch e
            row.note = sprint(showerror, e)[1:min(end, 80)]
            c = complex(NaN, NaN)
        end
        pts = sum(length, S.legs)
        converged
    end
    row.time = stats.time; row.bytes = stats.bytes
    row.iterations = it; row.points = pts; row.maxpoints = maxpts
    row.finish = stats.value ? "converged" : "maxiter"
    classify!(row, c, ref; tol_ok=tol_ok)
    if row.status == "ok" && !stats.value
        row.note *= " (maxiter)"
    end
    return row
end

"""
Current implementation (`Mandelbrot.spider_solve`), or with `mod=UpgradeSpider` the
shoulder-first version of 38caca9. Keyword arguments are passed on.
"""
function run_new(θ::Rational; impl::String="new", mod::Module=M, tol_ok=1e-10, kwargs...)
    ref = reference(θ)
    row = Row(impl=impl, angle=string(θ), preperiod=ref.preperiod, period=ref.period)
    local r
    GC.gc()
    stats = @timed try
        r = mod.spider_solve(θ; kwargs...)
        true
    catch e
        row.note = sprint(showerror, e)[1:min(end, 80)]
        false
    end
    row.time = stats.time; row.bytes = stats.bytes
    if !stats.value
        row.status = "exception"
        return row
    end
    row.iterations = r.iterations
    row.points = sum(length, r.spider.legs)
    row.finish = string(r.finish)
    row.note = "refine=$(r.stats.refinements) sat=$(r.stats.saturated) pruned=$(r.stats.pruned) newton=$(r.newton_attempts)"
    classify!(row, r.c, ref; tol_ok=tol_ok)
    if row.status == "ok" && !r.converged
        row.status = "nonconverged"
    end
    return row
end

end # module
