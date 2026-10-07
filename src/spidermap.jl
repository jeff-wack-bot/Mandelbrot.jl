#=
The spider algorithm of Hubbard and Schleicher.

We work with the normalisation used in [Hubbard_Schleicher_1995](@cite):

    p_λ(z) = λ(1 + z/2)^2,      c = λ/2,

so the critical point is z = -2, the critical value is z = 0 and p_λ(0) = λ.
A spider consists of `n` legs, one for each angle θ_i of the orbit of θ under
doubling. Leg `i` is a polyline from its shoulder (`legs[i][1]`, on the circle of
radius `SPIDER_RADIUS` at angle ≈ θ_i, continued radially to infinity) to its foot
(`legs[i][end]`). The feet are the marked points z_i; z_1 = 0 and z_2 = λ.

One step of the spider map replaces leg i by the lift under p_λ of leg σ(i), where
σ(i) is the index of 2θ_i in the orbit (`goestoidx`).

Implementation notes (see spiderlab/REPORT.md for the measurements behind them):

* Lifting a polyline by continuation is exact: a straight segment [a,b] that misses
  the origin subtends an angle < π at 0, so its continuous square root subtends an
  angle < π/2 and the correct branch at b is the one with Re(conj(√a)·√b) > 0. No
  branch cut is involved.
* Which of the two lifts becomes leg i is decided at infinity. Leg 1 is the lift of
  leg 2 whose foot is u = +1 (i.e. z = 0); every other leg is the lift whose
  shoulder direction differs from that of leg 1 by θ_i - θ_1. This is equivalent to
  the kneading-sequence rule but needs no sector bookkeeping.
* The chord between two lifted vertices need not be homotopic to the true lift of the
  segment. Segments are subdivided (before lifting) until no possible marked point
  lies between chord and arc.
* Vertices are removed when the triangle they span with their neighbours contains no
  marked point (this does not change the homotopy class of the leg), which keeps the
  number of points bounded.
* The feet can be carried in a different (higher) precision than the legs: the legs
  are only used to decide branches, the feet are lifted exactly.
=#

const SPIDER_RADIUS = 1.0e6

struct Spider{T<:AbstractFloat,L<:AbstractFloat}
    angle::Rational
    orbit::Sequence
    kneading_sequence::Sequence{KneadingSymbol{4}}
    legs::Vector{Vector{Complex{L}}}
    feet::Vector{Complex{T}}
end

function Base.show(io::IO, S::Spider{T,L}) where {T,L}
    npoints = sum(length, S.legs)
    print(io, "A ", repr(S.angle), "-spider with kneading sequence ", repr(S.kneading_sequence),
          " and ", npoints, " points (feet: ", T, ", legs: ", L, ")")
end

"""
    standardspider(theta; feet=Float64, legs=Float64)

The starting spider: radial legs at the angles of the orbit of `theta`, feet on the
unit circle, translated so that the foot of leg 1 sits at the critical value 0.
"""
function standardspider(theta::Rational; feet::Type{<:AbstractFloat}=Float64,
                        legs::Type{<:AbstractFloat}=Float64)
    orb = orbit(theta)
    K = thetaitinerary(theta, orb)
    T, L = feet, legs
    base = cispi(2 * T(orb.items[1]))
    F = [cispi(2 * T(a)) - base for a in orb.items]
    F[1] = zero(Complex{T})
    Ls = Vector{Vector{Complex{L}}}(undef, length(F))
    for (i, a) in enumerate(orb.items)
        shoulder = Complex{L}(SPIDER_RADIUS * cispi(2 * L(a)))
        Ls[i] = [shoulder, Complex{L}(F[i])]
    end
    return Spider{T,L}(theta, orb, K, Ls, F)
end

#-----------------------------------------------------------------------------
# Marked points with a cheap spatial query (sorted by real part).

struct PointSet{L}
    x::Vector{L}
    y::Vector{L}
end

function PointSet(pts::AbstractVector{Complex{L}}) where L
    p = sort(pts; by=real)
    return PointSet{L}(real.(p), imag.(p))
end

# smallest distance from any point of P to the segment [a,b], only looking at
# points within distance r of the segment's bounding box; returns Inf if none.
function mindist(P::PointSet{L}, a::Complex{L}, b::Complex{L}, r::L) where L
    xlo = min(real(a), real(b)) - r
    xhi = max(real(a), real(b)) + r
    ylo = min(imag(a), imag(b)) - r
    yhi = max(imag(a), imag(b)) + r
    best = L(Inf)
    i = searchsortedfirst(P.x, xlo)
    @inbounds while i <= length(P.x) && P.x[i] <= xhi
        y = P.y[i]
        if ylo <= y <= yhi
            q = complex(P.x[i], y)
            if q != a && q != b    # the segment's own end points are not obstacles
                best = min(best, segdist(q, a, b))
            end
        end
        i += 1
    end
    return best
end

function segdist(p::Complex, a::Complex, b::Complex)
    d = b - a
    l2 = abs2(d)
    iszero(l2) && return abs(p - a)
    t = clamp(real(conj(p - a) * d) / l2, 0, 1)
    return abs(p - (a + t * d))
end

cross2(a::Complex, b::Complex) = real(a) * imag(b) - imag(a) * real(b)

# does the closed triangle (a,b,c) contain a point of P other than its corners?
function anyintriangle(P::PointSet{L}, a::Complex{L}, b::Complex{L}, c::Complex{L}) where L
    xlo = min(real(a), real(b), real(c)); xhi = max(real(a), real(b), real(c))
    ylo = min(imag(a), imag(b), imag(c)); yhi = max(imag(a), imag(b), imag(c))
    o = cross2(b - a, c - a)
    i = searchsortedfirst(P.x, xlo)
    @inbounds while i <= length(P.x) && P.x[i] <= xhi
        y = P.y[i]
        if ylo <= y <= yhi
            p = complex(P.x[i], y)
            if p != a && p != b && p != c
                if iszero(o)
                    # degenerate triangle: only a point on the segment a–c blocks
                    segdist(p, a, c) == 0 && return true
                else
                    s1 = cross2(b - a, p - a) * o
                    s2 = cross2(c - b, p - b) * o
                    s3 = cross2(a - c, p - c) * o
                    (s1 >= 0 && s2 >= 0 && s3 >= 0) && return true
                end
            end
        end
        i += 1
    end
    return false
end

#-----------------------------------------------------------------------------
# Lifting

Base.@kwdef struct SpiderOptions
    refine::Bool = true        # subdivide segments whose lifted chord is not homotopic to the arc
    prune::Bool = true         # drop vertices whose removal does not change the homotopy class
    maxsweep::Float64 = 0.5    # max angle (rad) a segment may subtend at the critical value
    safety::Float64 = 2.0      # chord–arc deviation must be < dist(obstacle)/safety
    maxdepth::Int = 40         # recursion limit for subdivision of one segment
    maxpoints::Int = 10_000_000
end

mutable struct SpiderStats
    refinements::Int
    saturated::Int     # number of segments where maxdepth was hit
    pruned::Int
end
SpiderStats() = SpiderStats(0, 0, 0)

# continuous square root: the root of v nearest to (in the sense of the half plane
# around) the previous root u
@inline function nextroot(u::Complex, v::Complex)
    w = sqrt(v)
    return real(conj(u) * w) < 0 ? -w : w
end

@inline function sweep(va::Complex, vb::Complex)
    (iszero(va) || iszero(vb)) && return zero(real(va))
    return abs(angle(vb / va))
end

# Push the lift of the segment [va, vb] (u = √v, branch continuing from ua) onto Q.
function lift_segment!(Q, va, ua, vb, obst, opts::SpiderOptions, stats::SpiderStats, depth::Int)
    ub = nextroot(ua, vb)
    if opts.refine && depth < opts.maxdepth
        split = sweep(va, vb) > opts.maxsweep
        if !split
            um = nextroot(ua, (va + vb) / 2)
            dev = abs(um - (ua + ub) / 2)
            # a deviation below rounding cannot be resolved by splitting (this happens
            # when marked points merge, e.g. for Misiurewicz angles whose orbit angles
            # share a landing point)
            if dev > 4 * eps(real(typeof(ua))) * max(abs(ua), abs(ub))
                split = mindist(obst, ua, ub, opts.safety * dev) <= opts.safety * dev
            end
        end
        if split
            stats.refinements += 1
            vm = (va + vb) / 2
            lift_segment!(Q, va, ua, vm, obst, opts, stats, depth + 1)
            lift_segment!(Q, vm, Q[end], vb, obst, opts, stats, depth + 1)
            return Q
        end
    elseif opts.refine
        stats.saturated += 1
    end
    push!(Q, ub)
    return Q
end

# Lift the polyline P/λ under √, starting with the principal root at P[1].
function lift_path(P::Vector{Complex{L}}, λ::Complex{L}, obst, opts, stats) where L
    va = P[1] / λ
    ua = sqrt(va)
    Q = Complex{L}[ua]
    sizehint!(Q, length(P))
    for k in 2:length(P)
        vb = P[k] / λ
        lift_segment!(Q, va, ua, vb, obst, opts, stats, 0)
        va = vb
        ua = Q[end]
    end
    return Q
end

# Remove vertices whose triangle with the neighbours is free of marked points.
function prune_leg!(leg::Vector{Complex{L}}, marks::PointSet{L}, stats) where L
    n = length(leg)
    n <= 2 && return leg
    out = 1                  # leg[1:out] is the pruned prefix
    for k in 2:(n-1)
        a = leg[out]; b = leg[k]; c = leg[k+1]
        if anyintriangle(marks, a, b, c)
            out += 1
            leg[out] = b
        else
            stats.pruned += 1
        end
    end
    out += 1
    leg[out] = leg[n]
    resize!(leg, out)
    return leg
end

angdiff(a::Real, b::Real) = abs(mod(a - b + 0.5, 1) - 0.5)   # distance on R/Z

"""
    spider_map(S::Spider; options...)

One step of the spider map. See `SpiderOptions` for the keyword arguments.
"""
function spider_map(S::Spider{T,L}; stats::SpiderStats=SpiderStats(), kwargs...) where {T,L}
    opts = SpiderOptions(; kwargs...)
    n = length(S.legs)
    θs = S.orbit.items
    src = goestoidx(S.orbit)          # new leg i is a lift of old leg src[i]
    λ = S.feet[2]
    λL = Complex{L}(λ)

    # roots of z_k/λ in the u-plane: both preimages of every marked point.
    roots = [sqrt(Complex{L}(z) / λL) for z in S.feet]
    obst = PointSet(vcat(roots, -roots))

    # lift each needed source leg once (principal branch at the shoulder)
    lifts = Dict{Int,Vector{Complex{L}}}()
    for j in unique(src)
        lifts[j] = lift_path(S.legs[j], λL, obst, opts, stats)
    end

    # decide the branches
    signs = Vector{Int}(undef, n)
    U1 = lifts[src[1]]
    signs[1] = real(U1[end]) > 0 ? 1 : -1        # foot of leg 1 must be u = +1
    ψ1 = angle(signs[1] * U1[1]) / (2π)
    for i in 2:n
        U = lifts[src[i]]
        target = ψ1 + Float64(θs[i] - θs[1])
        ψ = angle(U[1]) / (2π)
        signs[i] = angdiff(ψ, target) < angdiff(ψ + 0.5, target) ? 1 : -1
    end

    # new feet, in the feet precision; branch from the leg's foot
    newfeet = Vector{Complex{T}}(undef, n)
    for i in 1:n
        z = S.feet[src[i]]
        if iszero(z)
            newfeet[i] = Complex{T}(-2)
        else
            r = sqrt(z / λ)
            ufoot = signs[i] * lifts[src[i]][end]
            if real(conj(Complex{L}(r)) * ufoot) < 0
                r = -r
            end
            newfeet[i] = 2 * (r - 1)
        end
    end
    newfeet[1] = zero(Complex{T})

    newlegs = Vector{Vector{Complex{L}}}(undef, n)
    for i in 1:n
        U = lifts[src[i]]
        s = signs[i]
        leg = Vector{Complex{L}}(undef, length(U) + 1)
        @inbounds for k in eachindex(U)
            leg[k+1] = 2 * (s * U[k] - 1)
        end
        leg[end] = Complex{L}(newfeet[i])
        # radial extension back out to the circle of radius SPIDER_RADIUS
        leg[1] = SPIDER_RADIUS * leg[2] / abs(leg[2])
        newlegs[i] = leg
    end

    if opts.prune
        marks = PointSet(Complex{L}.(newfeet))
        for leg in newlegs
            prune_leg!(leg, marks, stats)
        end
    end

    total = sum(length, newlegs)
    total > opts.maxpoints && error("spider has $total points (maxpoints = $(opts.maxpoints))")

    return Spider{T,L}(S.angle, S.orbit, S.kneading_sequence, newlegs, newfeet)
end

function spideriterates(S0::Spider, n_iter::Int; kwargs...)
    list = [S0]
    for i in 1:n_iter
        push!(list, spider_map(list[end]; kwargs...))
    end
    return list
end

#-----------------------------------------------------------------------------
# Newton polishing of the limit

# Nucleus: f_c^p(0) = 0.
function newton_nucleus(c::Complex{T}, p::Int; maxit::Int=100) where T
    tol = 64 * eps(T)
    for _ in 1:maxit
        z = zero(c); dz = zero(c)
        for _ in 1:p
            dz = 2z * dz + 1
            z = z^2 + c
        end
        iszero(z) && return c, true
        δ = z / dz
        c -= δ
        isfinite(c) || return c, false
        abs(δ) <= tol * max(1, abs(c)) && return c, true
    end
    return c, false
end

# Misiurewicz point: f^{q+p}(0) = f^q(0), deflated by the roots with smaller preperiod.
function newton_misiurewicz(c::Complex{T}, q::Int, p::Int; maxit::Int=100) where T
    tol = 64 * eps(T)
    zs = Vector{Complex{T}}(undef, q + p + 1)
    dzs = Vector{Complex{T}}(undef, q + p + 1)
    for _ in 1:maxit
        z = zero(c); dz = zero(c)
        zs[1] = z; dzs[1] = dz
        for k in 1:(q+p)
            dz = 2z * dz + 1
            z = z^2 + c
            zs[k+1] = z; dzs[k+1] = dz
        end
        g = zs[q+p+1] - zs[q+1]
        iszero(g) && return c, true
        r = (dzs[q+p+1] - dzs[q+1]) / g
        for i in 0:(q-1)
            h = zs[i+p+1] - zs[i+1]
            iszero(h) || (r -= (dzs[i+p+1] - dzs[i+1]) / h)
        end
        δ = 1 / r
        c -= δ
        isfinite(c) || return c, false
        abs(δ) <= tol * max(1, abs(c)) && return c, true
    end
    return c, false
end

function newton_polish(S::Spider, c)
    l = S.orbit.preperiod
    k = S.orbit.period
    return l == 0 ? newton_nucleus(c, k) : newton_misiurewicz(c, l + 1, k)
end

"""
    newton_feet(S::Spider; maxit=60)

Multiple-shooting Newton's method started from all feet of the spider: the unknowns are
c and the whole critical orbit Z_i = c(1 + z_i/2), the equations Z_{i+1} = Z_i^2 + c
together with the closing condition (Z_n = 0 for periodic angles, Z_n^2 + c = Z_{l+1}
for preperiodic ones). The bidiagonal structure is eliminated in O(n) per step. Starting
from the whole orbit instead of from c alone gives a much larger basin of attraction than
Newton's method for f_c^p(0) = 0. Returns (c, converged).
"""
function newton_feet(S::Spider{T}; maxit::Int=60) where T
    n = length(S.feet)
    l = S.orbit.preperiod
    c = S.feet[2] / 2
    iszero(c) && return c, false
    Z = [c * (1 + z / 2) for z in S.feet]
    Z[1] = c
    a = similar(Z); b = similar(Z)
    tol = 64 * eps(T)
    for _ in 1:maxit
        # δZ_i = a_i δc + b_i   (δZ_1 = δc)
        a[1] = one(c); b[1] = zero(c)
        for i in 1:(n-1)
            a[i+1] = 2Z[i] * a[i] + 1
            b[i+1] = 2Z[i] * b[i] + (Z[i]^2 + c - Z[i+1])
        end
        if l == 0
            # Z_n + δZ_n = 0
            δc = -(Z[n] + b[n]) / a[n]
        else
            # Z_{l+1} + δZ_{l+1} = (Z_n + δZ_n)^2 + c + δc
            δc = (Z[n]^2 + c - Z[l+1] - b[l+1] + 2Z[n] * b[n]) / (a[l+1] - 2Z[n] * a[n] - 1)
        end
        isfinite(δc) || return c, false
        for i in 1:n
            Z[i] += a[i] * δc + b[i]
        end
        c = Z[1]
        abs(δc) <= tol * max(1, abs(c)) && return c, true
    end
    return c, false
end

"""
    component_size(S::Spider, c)

Length scale on which other solutions of the same equation can be expected near c:
the nucleus size estimate 1/|β λ²| for periodic angles, and 1/|∂_c f_c^q(0)| for
preperiodic ones (q = preperiod of the critical point).
"""
function component_size(S::Spider, c)
    l = S.orbit.preperiod
    k = S.orbit.period
    if l == 0
        z = c; λ = one(c); β = one(c)
        for _ in 1:(k-1)
            λ *= 2z
            β += 1 / λ
            z = z^2 + c
        end
        return Float64(1 / abs(β * λ^2))
    else
        z = zero(c); dz = zero(c)
        for _ in 1:(l+1)
            dz = 2z * dz + 1
            z = z^2 + c
        end
        return Float64(1 / abs(dz))
    end
end

"""
    exact_type(S::Spider, c)

Whether `c` has exactly the combinatorial type of the spider's angle: for a periodic
angle of period p, f_c^d(0) ≠ 0 for d < p (Newton's method also finds nuclei whose
period divides p, e.g. c = 0); for a preperiodic angle (preperiod l, period k) the
critical value is not periodic and f_c^l(c) is not already on the cycle.
"""
function exact_type(S::Spider, c)
    isfinite(c) && !iszero(c) || return false
    n = length(S.feet)
    l = S.orbit.preperiod
    k = S.orbit.period
    Z = Vector{typeof(c)}(undef, n + 1)
    Z[1] = c
    for i in 1:n
        Z[i+1] = Z[i]^2 + c
    end
    tiny = 1e3 * eps(real(typeof(c))) * max(1, maximum(abs, Z))
    if l == 0
        return all(abs(Z[d]) > tiny for d in 1:(k-1))
    else
        all(abs(Z[d]) > tiny for d in 1:n) || return false
        return l == 0 || abs(Z[l] - Z[l+k]) > tiny
    end
end

# The marked points of the spider belonging to the parameter c: z_i = 2(f_c^{i-1}(c)/c - 1).
function critical_feet(c::Complex{T}, n::Int) where T
    F = Vector{Complex{T}}(undef, n)
    Z = c
    for i in 1:n
        F[i] = 2 * (Z / c - 1)
        Z = Z^2 + c
    end
    F[1] = zero(Complex{T})
    return F
end

# do the closed segments [p1,p2] and [q1,q2] intersect?
function segments_intersect(p1::Complex, p2::Complex, q1::Complex, q2::Complex)
    d1 = cross2(p2 - p1, q1 - p1)
    d2 = cross2(p2 - p1, q2 - p1)
    d3 = cross2(q2 - q1, p1 - q1)
    d4 = cross2(q2 - q1, p2 - q1)
    if ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))
        return true
    end
    # collinear / touching cases
    onseg(a, b, x) = min(real(a), real(b)) <= real(x) <= max(real(a), real(b)) &&
                     min(imag(a), imag(b)) <= imag(x) <= max(imag(a), imag(b))
    (d1 == 0 && onseg(p1, p2, q1)) && return true
    (d2 == 0 && onseg(p1, p2, q2)) && return true
    (d3 == 0 && onseg(q1, q2, p1)) && return true
    (d4 == 0 && onseg(q1, q2, p2)) && return true
    return false
end

# Winding numbers, around each point of F except F[i], of the closed curve: leg A from
# its foot back to its shoulder, the shorter arc at infinity to the shoulder of leg B,
# leg B back to the (common) foot.
function leg_windings(A::Vector{<:Complex}, B::Vector{<:Complex}, F::Vector{<:Complex}, i::Int)
    darc = mod(angle(B[1]) - angle(A[1]) + π, 2π) - π
    w = zeros(Int, length(F))
    for j in eachindex(F)
        j == i && continue
        z = F[j]
        tot = darc
        for k in length(A):-1:2          # A reversed: foot -> shoulder
            tot += angle((A[k-1] - z) / (A[k] - z))
        end
        for k in 1:(length(B)-1)         # B: shoulder -> foot
            tot += angle((B[k+1] - z) / (B[k] - z))
        end
        isfinite(tot) || return nothing   # a leg runs through a foot (unresolvable in this precision)
        w[j] = round(Int, tot / 2π)
    end
    return w
end

"""
    certify(S::Spider, c; drag=:isotopy, legs=true, steps=1, kwargs...)

Check that `c` is the parameter the spider `S` is converging to. The feet of `S` are
moved to the marked points Z_i = 2(f_c^{i-1}(c)/c - 1) of `c`, giving a spider T. With
`drag=:isotopy` the feet travel along straight paths that must not meet other legs or
other paths (so T is obtained by an isotopy); with `drag=:room` every foot must move
less than a quarter of its distance to the other feet and legs. Then the spider map is
applied `steps` times. `c` is accepted if

* for every foot the legs select the preimage that is the actual orbit point of `c`
  (so the feet of T are fixed), and
* with `legs=true`: every lifted leg is homotopic to the old one, tested by the winding
  numbers of (old leg ∪ arc at ∞ ∪ new leg) around the other feet (all equal, since arcs
  to ∞ are only defined up to a full twist).

Together these say that T is (numerically) a fixed point of the spider map, which by
Thurston rigidity is unique. Returns `:certified`, `:rejected` or `:inconclusive`
(paths cross legs, or feet have merged as happens for Misiurewicz angles whose orbit
angles share a landing point).
"""
function certify(S::Spider{T,L}, c; drag::Symbol=:isotopy, legs::Bool=true, steps::Int=1,
                 kwargs...) where {T,L}
    n = length(S.feet)
    Z = critical_feet(Complex{T}(c), n)
    scale = max(1, maximum(abs, Z))
    ZL = Complex{L}.(Z)
    FL = Complex{L}.(S.feet)
    for i in 1:n, j in (i+1):n
        abs(Z[i] - Z[j]) <= 1e3 * eps(T) * scale && return :inconclusive   # merged feet
    end
    if drag === :isotopy
        for i in 1:n
            a, b = FL[i], ZL[i]
            a == b && continue
            for j in 1:n
                j == i && continue
                leg = S.legs[j]
                for k in 1:(length(leg)-1)
                    segments_intersect(a, b, leg[k], leg[k+1]) && return :inconclusive
                end
                if j > i && FL[j] != ZL[j]
                    segments_intersect(a, b, FL[j], ZL[j]) && return :inconclusive
                end
            end
        end
    else
        for i in 1:n
            room = Inf
            for j in 1:n
                j == i && continue
                room = min(room, Float64(abs(FL[j] - FL[i])))
                leg = S.legs[j]
                for k in 1:(length(leg)-1)
                    room = min(room, Float64(segdist(FL[i], leg[k], leg[k+1])))
                end
            end
            Float64(abs(ZL[i] - FL[i])) < room / 4 || return :inconclusive
        end
    end
    S2 = Spider{T,L}(S.angle, S.orbit, S.kneading_sequence,
                     [push!(copy(leg), ZL[i]) for (i, leg) in enumerate(S.legs)], Z)
    for _ in 1:steps
        S3 = spider_map(S2; kwargs...)
        # Z[σ(i)] has the preimages Z[i] and -Z[i]-4; the legs must select Z[i]
        # (for the critical point -2 both preimages coincide, so the choice is immaterial)
        ok = all(abs(S3.feet[i] - Z[i]) <= abs(-S3.feet[i] - 4 - Z[i]) ||
                 abs(S3.feet[i] + 2) <= 1e3 * eps(T) * scale for i in 1:n)
        ok || return :rejected
        newlegs = S3.legs
        for i in 1:n
            newlegs[i][end] = ZL[i]
        end
        if legs
            for i in 1:n
                w = leg_windings(S2.legs[i], newlegs[i], ZL, i)
                w === nothing && return :inconclusive
                ws = [w[j] for j in 1:n if j != i]
                isempty(ws) || allequal(ws) || return :rejected
            end
        end
        S2 = Spider{T,L}(S3.angle, S3.orbit, S3.kneading_sequence, newlegs, Z)
    end
    return :certified
end

#-----------------------------------------------------------------------------

"""
    SpiderResult

Returned by `spider_solve`: `c` the parameter, `converged`, `iterations` of the
spider map, how the last digits were obtained (`finish`: `:spider`, `:newton_certified`
or `:newton_trusted`), the final spider, the step history |c_k - c_{k-1}| and some
counters.
"""
struct SpiderResult{C,S}
    c::C
    converged::Bool
    iterations::Int
    finish::Symbol
    spider::S
    history::Vector{Float64}
    stats::SpiderStats
    newton_attempts::Int
    size::Float64      # see `component_size`; if size < eps(T)|c| the precision cannot resolve the answer
end

"""
    spider_solve(theta; maxiter=2000, tol=8eps(feet), feet=Float64, legs=Float64,
                 polish=true, try_newton_at=Inf, trust_fraction=1e-3, kwargs...)

Run the spider algorithm for the angle `theta` and return a [`SpiderResult`](@ref).
The feet (landing points) are carried in precision `feet`, the legs in precision
`legs`. `feet=BigFloat, legs=Float64` is the recommended way to get parameters beyond
Float64 accuracy: it costs about 1.3× a Float64 run, while BigFloat legs cost ~150×.

With `polish=true`, whenever the step |c_k - c_{k-1}| has halved since the last attempt
(and is below `try_newton_at`, relative), a multiple-shooting Newton step
([`newton_feet`](@ref)) is started from the current feet. Its root is accepted if it has
the exact combinatorial type of `theta` and [`certify`](@ref) confirms it. If
certification is inconclusive, the root is accepted only when the spider is visibly
converging to it: the last steps decrease geometrically and the a posteriori estimate
Δ r/(1-r) is below `trust_fraction` times [`component_size`](@ref). Otherwise the plain
spider iteration runs until the step is below `tol`.

Check `result.size` against `eps(feet)*abs(result.c)`: if the component is smaller than
the working precision, the answer is accurate to that precision but does not identify
the component (use `feet=BigFloat`).

Remaining keyword arguments are passed to `spider_map` (see `SpiderOptions`).
"""
function spider_solve(theta::Rational; feet::Type{<:AbstractFloat}=Float64,
                      legs::Type{<:AbstractFloat}=Float64, kwargs...)
    θs = orbit(theta).items
    gap = minimum((abs(mod(a - b + 1//2, 1) - 1//2) for a in θs for b in θs if a != b); init=1//2)
    if gap < 1e3 * eps(legs)
        @warn "orbit angles of $theta are only $(Float64(gap)) apart; legs in $legs cannot separate their shoulders. Use legs=BigFloat (slow) for reliable results."
    end
    return spider_solve(standardspider(theta; feet=feet, legs=legs); kwargs...)
end

function spider_solve(S::Spider{T,L}; maxiter::Int=2000, tol::Real=8 * eps(T), polish::Bool=true,
                      try_newton_at::Real=Inf, trust_fraction::Real=1e-3, kwargs...) where {T,L}
    stats = SpiderStats()
    history = Float64[]
    attempts = 0
    c_last = S.feet[2] / 2
    if S.orbit.items == [0//1]   # θ = 0: the main cardioid, c = 0
        return SpiderResult(zero(c_last), true, 0, :spider, S, history, stats, 0, Inf)
    end
    Δ_try = Inf
    size = NaN
    for it in 1:maxiter
        S = spider_map(S; stats=stats, kwargs...)
        c = S.feet[2] / 2
        Δ = Float64(abs(c - c_last))
        push!(history, Δ)
        isfinite(Δ) || return SpiderResult(c, false, it, :spider, S, history, stats, attempts, size)
        scale = max(1.0, Float64(abs(c)))
        # try Newton whenever the spider has moved substantially closer since the last try
        if polish && it >= 2 && Δ <= try_newton_at * scale && Δ <= Δ_try / 2
            Δ_try = Δ
            attempts += 1
            cp, ok = newton_feet(S)
            if ok && exact_type(S, cp)
                verdict = certify(S, cp; kwargs...)
                if verdict === :certified
                    return SpiderResult(cp, true, it, :newton_certified, S, history, stats, attempts,
                                        component_size(S, cp))
                elseif verdict === :inconclusive
                    # feet too close together to certify: trust Newton only once the spider
                    # itself is well inside the basin of this solution
                    # (one small step proves nothing: require a steady geometric decrease
                    # and use the a posteriori estimate Δ r/(1-r) of the distance to the limit)
                    size = component_size(S, cp)
                    m = length(history)
                    steady = m >= 5 && all(history[k] < history[k-1] for k in (m-3):m)
                    r = steady ? history[m] / history[m-1] : 1.0
                    est = r < 0.99 ? Δ * r / (1 - r) : Inf
                    if isfinite(size) && est <= trust_fraction * size && abs(cp - c) <= 2est + 100 * eps(T) * scale
                        return SpiderResult(cp, true, it, :newton_trusted, S, history, stats, attempts, size)
                    end
                    Δ_try = 2Δ   # retry at the next step
                end
            end
        end
        if Δ <= tol * scale
            return SpiderResult(c, true, it, :spider, S, history, stats, attempts, component_size(S, c))
        end
        c_last = c
    end
    c = S.feet[2] / 2
    return SpiderResult(c, false, maxiter, :spider, S, history, stats, attempts, component_size(S, c))
end

"""
    parameter(S0::Spider, max_iter::Int; kwargs...)
    parameter(theta::Rational, max_iter::Int=1000; kwargs...)

Computes the parameter (center) of a hyperbolic component using the spider algorithm
described in [Hubbard_Schleicher_1995](@cite). For a preperiodic angle the result is
the corresponding Misiurewicz parameter.

Keyword arguments are passed to [`spider_solve`](@ref); for example
`parameter(1//7; feet=BigFloat)` returns the parameter in extended precision.

# References
- [Hubbard_Schleicher_1995](@cite): The spider algorithm
"""
function parameter(S::Spider, max_iter::Int; kwargs...)
    return spider_solve(S; maxiter=max_iter, kwargs...).c
end

function parameter(theta::Rational, max_iter::Int=1000; kwargs...)
    return spider_solve(theta; maxiter=max_iter, kwargs...).c
end
