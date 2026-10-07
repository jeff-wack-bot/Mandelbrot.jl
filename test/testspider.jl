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
