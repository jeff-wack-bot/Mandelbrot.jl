# Spider Algorithm

The spider algorithm computes the parameter of a quadratic polynomial from an external angle [Hubbard_Schleicher_1995](@cite). It works by iterating a "spider" — a collection of legs connecting the critical point to infinity — until convergence.

## ``2/9``

```@example spider
using Mandelbrot
using BrotViz
using CairoMakie

angle = 2//9
S0 = Mandelbrot.standardspider(angle)
list = Mandelbrot.spideriterates(S0, 20)

fig = Figure()
ax = Axis(fig[1, 1])

record(fig, "spider_2_9.gif", 1:20; framerate=3) do i
    empty!(ax)
    spiderplot!(ax, list[i])
end
nothing # hide
```

![Spider algorithm for 2/9](spider_2_9.gif)

## Introduction

```@example standardspider
using Mandelbrot
using BrotViz
using CairoMakie

theta = 9//56
S = Mandelbrot.standardspider(theta)

fig = Figure()
ax = Axis(fig[1,1])

spiderplot!(ax,S)

fig
```

```@example kneadingseq
using Mandelbrot

K = KneadingSequence(4//15)
```





## Accuracy, precision and options

`parameter(θ)` runs [`Mandelbrot.spider_solve`](@ref), which returns more information:

```@example spidersolve
using Mandelbrot
r = Mandelbrot.spider_solve(1//1023)
(r.c, r.iterations, r.finish, r.size)
```

* The legs are kept minimal: every vertex whose removal does not change the homotopy
  class of the leg (relative to the marked points) is pruned, and segments are
  subdivided before lifting whenever the lifted chord could pass on the wrong side of a
  marked point. A spider therefore has a few points per leg instead of growing every
  iteration.
* Once the spider is roughly in place, Newton's method on the whole critical orbit
  (multiple shooting) proposes the answer, and [`Mandelbrot.certify`](@ref) checks it
  against the spider. Slowly converging angles (near the cusp, `1//(2^p-1)`) finish in a
  handful of iterations instead of thousands.
* `feet=BigFloat` carries the landing points in extended precision while the legs stay
  in `Float64`. This is needed when the component is smaller than Float64 can resolve
  (compare `r.size` with `eps()*abs(r.c)`), e.g. real components next to ``c = -2``:

```@example spidersolve
θ = (big(2)^29 - 1) // (big(2)^30 - 1)
r = Mandelbrot.spider_solve(θ; feet=BigFloat)
(Float64(real(r.c)), r.size)
```

For a pure `Float64` run the same angle returns a number within `1e-16` of the center,
but the component (size ``\sim 10^{-34}``) cannot be identified at that precision.

```@docs
Mandelbrot.spider_solve
Mandelbrot.certify
```

## References

```@bibliography
Pages = [@__FILE__]
```

