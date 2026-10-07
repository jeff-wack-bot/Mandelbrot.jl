# Summaries of results/*.csv
include(joinpath(@__DIR__, "SpiderBench.jl")); using .SpiderBench
using Statistics, Printf
const STATUSES = ["ok", "inaccurate", "unresolved", "wrong", "nan", "nonconverged", "exception"]
function summarize(rows, key; impls=unique(r.impl for r in rows))
    keys_ = sort(unique(key.(rows)))
    @printf("%-14s %-10s %5s  %s  %9s %9s %10s %9s\n", "impl", "group", "n", join(rpad.(first.(STATUSES, 5), 6)), "med.it", "med.ms", "med.MB", "med.pts")
    for impl in impls, k in keys_
        rs = filter(r -> r.impl == impl && key(r) == k, rows)
        isempty(rs) && continue
        cnt = [count(r -> r.status == s, rs) for s in STATUSES]
        @printf("%-14s %-10s %5d  %s  %9.0f %9.2f %10.2f %9.0f\n", impl, string(k), length(rs), join(rpad.(string.(cnt), 6)),
                median(r.iterations for r in rs), 1e3 * median(r.time for r in rs), median(r.bytes for r in rs) / 2^20,
                median(r.points for r in rs))
    end
end
