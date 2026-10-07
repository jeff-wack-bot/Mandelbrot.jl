# spiderlab — testing, benchmarking and profiling the spider algorithm

A separate Julia environment (it `dev`s the package in `..`) used to produce
[REPORT.md](REPORT.md).

```
julia --project=spiderlab -e 'using Pkg; Pkg.instantiate()'
```

| file | purpose |
|---|---|
| `LegacySpider.jl` | verbatim copy of the spider implementation at upstream `619df42`, wrapped in a module. `LegacySpider.FIX_GETINDEX[] = true` switches on the corrected `Sequence` indexing. |
| `GroundTruth.jl` | independent reference values: parameter-ray tracing (`exray_in`) + Newton for the nucleus / Misiurewicz point, all in BigFloat. Shares no code with the spider. |
| `SpiderBench.jl` | harness: angle families (`periodic_angles`, `preperiodic_angles`, `family`), runners `run_legacy` / `run_new` (time, bytes, points, iterations), classification against the ground truth, CSV I/O. |
| `analyze.jl` | `summarize(rows, key)` tables. |
| `figures.jl` | renders `figures/*.png` from `results/*.csv`. |
| `experiments/sweep_small.jl` | every periodic angle with period ≤ 10 and every preperiodic angle with denominator ≤ 64, all implementations. |
| `experiments/scaling.jl` | families × period up to 64: legacy, Float64, BigFloat feet, BigFloat feet+legs. |
| `experiments/history_ablation.jl` | error / point-count histories; refinement and pruning ablation. |
| `experiments/certify_soundness*.jl` | adversarial test of `certify`: every wrong nucleus / Misiurewicz point offered as a candidate. |
| `experiments/profile.jl` | BenchmarkTools + Profile of `spider_map`. |

Status classes used everywhere: **ok** (|c − c_ref| ≤ 1e-10·max(1,|c|) and inside the right
component), **inaccurate** (right component, error above 1e-10), **unresolved** (as accurate
as Float64 allows, but the component is smaller than that), **wrong** (another component /
point), **nan**, **exception**.

All runs in REPORT.md were done on a Raspberry Pi 5 (4 × Cortex-A76, 16 GB), Julia 1.13.1,
single-threaded.
