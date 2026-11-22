using Mandelbrot, Test

# ## The spider algorithm
#
# `parameter(θ)` returns the center of the hyperbolic component whose root has
# external angle θ (θ periodic), or the Misiurewicz point where the ray of angle θ
# lands (θ preperiodic). Reference values below come from tracing the parameter ray
# and polishing with Newton's method in BigFloat (see `spiderlab/GroundTruth.jl`).

const REFERENCE = [
    (1//3,   -1.0 + 0.0im),
    (1//7,   -0.12256116687665362 + 0.7448617666197442im),   # rabbit
    (2//7,   -0.12256116687665362 + 0.7448617666197442im),
    (3//7,   -1.7548776662466927 + 0.0im),                   # airplane
    (1//15,   0.2822713907669139 + 0.5300606175785253im),
    (2//9,   -0.16359826155202264 + 1.0977806428882722im),
    (1//31,   0.3795135880159237 + 0.3349323055974976im),
    (11//31, -1.2563679300681807 + 0.3803209634727225im),
    (1//127,  0.37600868184676756 + 0.14474937132163287im),
    (1//6,    0.0 + 1.0im),                                  # Misiurewicz points
    (1//4,   -0.22815549365396182 + 1.1151425080399373im),
    (9//56,  -0.10109636384562216 + 0.9562865108091415im),
    (5//12,  -1.5436890126920764 + 0.0im),
    (1//2,   -2.0 + 0.0im),
]

for (θ, c) in REFERENCE
    @test isapprox(parameter(θ), c; atol=1e-13)
end

# Preperiodic kneading sequences index into their periodic part correctly:

K = KneadingSequence(1//6)
@test string.([K[i] for i in 1:5]) == ["A", "A", "B", "A", "B"]

# The spider stays small: pruning removes every vertex that does not change the
# homotopy class of a leg.

S = Mandelbrot.standardspider(1//1023)
for _ in 1:50
    global S = Mandelbrot.spider_map(S)
end
@test sum(length, S.legs) < 100

# The lift of the segment ``[1+i, 1-i]`` under ``\sqrt{\cdot}`` is an arc of the hyperbola
# ``x^2 - y^2 = 1`` through ``u = 1``, and its chord lies at ``x = \mathrm{Re}\sqrt{1+i}
# \approx 1.099``. A foot between the two is enclosed, one on either side is not:

va, vb = 1.0 + 1.0im, 1.0 - 1.0im
ua, ub = sqrt(va), sqrt(vb)
@test Mandelbrot.inlens(1.05 + 0.0im, ua, ub, va, vb)
@test !Mandelbrot.inlens(1.2 + 0.0im, ua, ub, va, vb)
@test !Mandelbrot.inlens(0.95 + 0.0im, ua, ub, va, vb)

# Each step lifts the legs into the right homotopy class: it agrees with the lift of the
# same spider whose legs have been subdivided into short pieces, so that every chord is
# homotopic to its arc.

function densify!(out, a, b)    # pieces shorter than 1e-3 of their distance to 0
    if abs(b - a) > 1e-3 * min(abs(a), abs(b)) + 1e-9
        m = (a + b) / 2
        densify!(out, a, m)
        densify!(out, m, b)
    else
        push!(out, b)
    end
    return out
end
densify(P) = foldl((out, k) -> densify!(out, P[k], P[k+1]), 1:length(P)-1; init=[P[1]])
S = Mandelbrot.standardspider(1//1023)
for _ in 1:12
    D = Mandelbrot.Spider{Float64,Float64}(S.angle, S.orbit, S.kneading_sequence,
                                           densify.(S.legs), copy(S.feet))
    R = Mandelbrot.spider_map(D; refine=false, prune=false)
    global S = Mandelbrot.spider_map(S)
    @test S.feet ≈ R.feet
    for i in eachindex(S.feet)
        w = Mandelbrot.leg_windings(S.legs[i], R.legs[i], S.feet, i)
        @test allequal(w[j] for j in eachindex(w) if j != i)
    end
end

# A real parameter with extended-precision feet: the Newton candidate lies within 1e-170 of
# the real axis, and the certificate must still go through.

r = Mandelbrot.spider_solve(3//5; feet=BigFloat)
@test r.finish === :newton_certified
@test isapprox(r.c, -1.3107026413368328; atol=1e-14)

# A slowly converging angle (near the cusp) is finished by a certified Newton step:

r = Mandelbrot.spider_solve(1//(2^20 - 1))
@test r.converged && r.finish === :newton_certified

# The feet (landing points) can be carried in extended precision while the legs stay
# in Float64. Here the component has diameter ~1e-29, far below Float64 resolution:

θ = (big(2)^25 - 1) // (big(2)^26 - 1)
r = Mandelbrot.spider_solve(θ; feet=BigFloat)
@test r.converged
@test r.size < 1e-25
z = zero(r.c)
for _ in 1:26
    global z = z^2 + r.c
end
@test abs(z) < 1e-60    # r.c is a nucleus of period 26 to high accuracy
