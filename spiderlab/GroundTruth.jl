# Independent reference values for the spider algorithm.
#
# For an external angle θ we trace the parameter ray R_M(θ) inwards with Newton's
# method (the classic "exray_in" scheme), then polish the end point with Newton's
# method for either the nucleus (θ periodic) or the Misiurewicz point (θ
# preperiodic). Everything runs in BigFloat. Nothing here shares code with the
# spider implementation.
module GroundTruth

export reference_parameter, exray_in, nucleus, misiurewicz, angle_type, size_estimate

"""Return (preperiod, period) of θ under doubling (exact rational arithmetic)."""
function angle_type(θ::Rational)
    θ = Rational{BigInt}(θ)
    seen = Dict{Rational{BigInt},Int}()
    k = 0
    while !haskey(seen, θ)
        seen[θ] = k
        θ = mod(2θ, 1)
        k += 1
    end
    return seen[θ], k - seen[θ]
end

"""
    exray_in(θ, depth; sharpness=8, er=2^16, prec=256)

Trace the parameter ray of angle θ from potential log(er) down to depth `depth`
(number of doublings). Returns the final point (Complex{BigFloat}).
"""
function exray_in(θ::Rational, depth::Int; sharpness::Int=8, er=big(2)^16, prec::Int=256)
    setprecision(BigFloat, prec) do
        φ = Rational{BigInt}(θ)
        R = BigFloat(er)
        c = R * cispi(2 * BigFloat(φ))
        for m in 1:depth
            for j in 1:sharpness
                r = R^(BigFloat(1//2)^(BigFloat(j) / sharpness))
                t = r * cispi(2 * BigFloat(φ))
                for _ in 1:64
                    z = zero(c); dz = zero(c)
                    for _ in 1:m
                        dz = 2z * dz + 1
                        z = z^2 + c
                    end
                    δ = (z - t) / dz
                    c -= δ
                    abs(δ) < eps(BigFloat) * 2^20 * max(1, abs(c)) && break
                end
            end
            φ = mod(2φ, 1)
        end
        return c
    end
end

"""Newton's method for f_c^p(0) = 0. Returns (c, converged::Bool)."""
function nucleus(c0, p::Int; prec::Int=256, maxit::Int=200)
    setprecision(BigFloat, prec) do
        c = Complex{BigFloat}(c0)
        for _ in 1:maxit
            z = zero(c); dz = zero(c)
            for _ in 1:p
                dz = 2z * dz + 1
                z = z^2 + c
            end
            iszero(z) && return c, true
            δ = z / dz
            c -= δ
            if abs(δ) <= eps(BigFloat) * 2^8 * max(1, abs(c))
                return c, true
            end
        end
        return c, false
    end
end

"""
Newton's method for the Misiurewicz equation f^{q+p}(0) = f^q(0), deflated by the
roots with smaller preperiod. q = preperiod of the critical point (= angle preperiod + 1).
"""
function misiurewicz(c0, q::Int, p::Int; prec::Int=256, maxit::Int=300)
    setprecision(BigFloat, prec) do
        c = Complex{BigFloat}(c0)
        for _ in 1:maxit
            zs = Vector{Complex{BigFloat}}(undef, q + p + 1)
            dzs = similar(zs)
            z = zero(c); dz = zero(c)
            zs[1] = z; dzs[1] = dz
            for k in 1:(q+p)
                dz = 2z * dz + 1
                z = z^2 + c
                zs[k+1] = z; dzs[k+1] = dz
            end
            g = zs[q+p+1] - zs[q+1]
            iszero(g) && return c, true
            dg = dzs[q+p+1] - dzs[q+1]
            r = dg / g
            for i in 0:(q-1)   # deflate f^{i+p} = f^i for i < q
                h = zs[i+p+1] - zs[i+1]
                dh = dzs[i+p+1] - dzs[i+1]
                iszero(h) && continue
                r -= dh / h
            end
            δ = 1 / r
            c -= δ
            if abs(δ) <= eps(BigFloat) * 2^8 * max(1, abs(c))
                return c, true
            end
        end
        return c, false
    end
end

"""Exact period check: smallest k with |f^k(0)| tiny (relative)."""
function exact_period(c, pmax::Int; tol=nothing)
    tol === nothing && (tol = sqrt(eps(real(typeof(c)))))
    z = zero(c)
    for k in 1:pmax
        z = z^2 + c
        abs(z) < tol && return k
    end
    return 0
end

"""
Nucleus size estimate (Heiland-Allen): for a period p nucleus c,
size ≈ 1 / |β λ²| with λ = ∏ 2z_k, β = Σ 1/∏_{j≤k} 2z_j  (k = 1..p-1).
"""
function size_estimate(c, p::Int)
    z = c; l = one(c); b = one(c)
    for _ in 1:(p-1)
        l *= 2z
        b += 1 / l
        z = z^2 + c
    end
    return 1 / abs(b * l^2)
end

"""
    reference_parameter(θ; prec=256, depthfactor=4, sharpness=8)

Returns a NamedTuple (c, ok, kind, preperiod, period, info).
"""
function reference_parameter(θ::Rational; prec::Int=256, depthfactor::Int=4, sharpness::Int=8,
                             mindepth::Int=24)
    l, p = angle_type(θ)
    if l == 0
        depth = max(mindepth, depthfactor * p)
        cr = exray_in(θ, depth; sharpness=sharpness, prec=prec)
        c, ok = nucleus(cr, p; prec=prec)
        per = setprecision(() -> exact_period(c, p; tol=BigFloat(2)^(-prec ÷ 2)), BigFloat, prec)
        ok &= per == p
        return (c=c, ok=ok, kind=:periodic, preperiod=0, period=p, ray=cr)
    else
        depth = max(mindepth, depthfactor * (l + p))
        cr = exray_in(θ, depth; sharpness=sharpness, prec=prec)
        c, ok = misiurewicz(cr, l + 1, p; prec=prec)
        return (c=c, ok=ok, kind=:preperiodic, preperiod=l, period=p, ray=cr)
    end
end

end # module
