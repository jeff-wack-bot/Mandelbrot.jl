# Figures for REPORT.md.  julia --project=spiderlab spiderlab/figures.jl
include(joinpath(@__DIR__, "SpiderBench.jl")); using .SpiderBench
using CairoMakie, Statistics

const FIG = joinpath(@__DIR__, "figures"); mkpath(FIG)
# reference palette (dataviz skill): categorical slots, status colours, sequential blue
const SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
const STATUS = Dict("ok" => "#0ca30c", "inaccurate" => "#fab219", "unresolved" => "#ec835a",
                    "wrong" => "#d03b3b", "nan" => "#d03b3b", "nonconverged" => "#ec835a",
                    "exception" => "#d03b3b")
const MARK = Dict("ok" => :circle, "inaccurate" => :diamond, "unresolved" => :utriangle,
                  "wrong" => :xcross, "nan" => :cross, "nonconverged" => :dtriangle, "exception" => :cross)
const INK = "#52514e"
set_theme!(Theme(fontsize=13, textcolor="#0b0b0b",
                 Axis=(xgridcolor=(:black, 0.06), ygridcolor=(:black, 0.06), spinewidth=0.8,
                       leftspinecolor=INK, bottomspinecolor=INK, topspinevisible=false,
                       rightspinevisible=false)))

function mset!(ax)
    xs = range(-2.1, 0.6, length=700); ys = range(-1.25, 1.25, length=650)
    img = [begin
        c = complex(x, y); z = 0im; n = 0
        while abs2(z) < 4 && n < 200; z = z^2 + c; n += 1; end
        n == 200 ? 1.0 : 0.0
    end for x in xs, y in ys]
    heatmap!(ax, xs, ys, img; colormap=[:white, "#e9e8e4"])
end

refc(r) = complex(r.c_re, r.c_im)
function refpoint(r)
    c = ComplexF64(reference(eval(Meta.parse(r.angle))).c)
    return Point2f(real(c), imag(c))
end

function fig_map()
    rows = read_csv(joinpath(@__DIR__, "results", "sweep_small.csv"))
    fig = Figure(size=(1250, 500))
    for (k, (impl, title)) in enumerate((("legacy", "legacy (as released)"), ("legacy+fix", "legacy + Sequence indexing fix"),
                                         ("new", "upgraded")))
        ax = Axis(fig[1, k], title=title, aspect=DataAspect(), xlabel="Re c", ylabel=k == 1 ? "Im c" : "")
        mset!(ax)
        rs = filter(r -> r.impl == impl, rows)
        for st in ("ok", "inaccurate", "wrong", "nan")
            sel = filter(r -> r.status == st, rs)
            isempty(sel) && continue
            pts = refpoint.(sel)
            pts = vcat(pts, [Point2f(p[1], -p[2]) for p in pts])   # θ and 1-θ
            scatter!(ax, pts; color=STATUS[st], marker=MARK[st], markersize=st == "ok" ? 5 : 9,
                     strokewidth=st == "ok" ? 0 : 0.5, strokecolor=:white, label="$st ($(length(sel)))")
        end
        axislegend(ax; position=:lb, framevisible=false, labelsize=11)
    end
    save(joinpath(FIG, "map_status.png"), fig, px_per_unit=2)

    # error of the legacy implementation for periodic angles
    fig = Figure(size=(900, 520))
    ax = Axis(fig[1, 1], title="legacy, periodic angles p ≤ 10: log₁₀ |c − c_ref| after ≤ 250 iterations",
              aspect=DataAspect(), xlabel="Re c", ylabel="Im c")
    mset!(ax)
    rs = filter(r -> r.impl == "legacy" && r.preperiod == 0, rows)
    pts = refpoint.(rs); e = log10.(max.([r.err for r in rs], 1e-17))
    pts2 = vcat(pts, [Point2f(p[1], -p[2]) for p in pts]); e2 = vcat(e, e)
    o = sortperm(e2)
    sc = scatter!(ax, pts2[o]; color=e2[o], colormap=cgrad(["#cde2fb", "#6da7ec", "#2a78d6", "#184f95", "#0d366b"]),
                  colorrange=(-16, -3), markersize=6)
    Colorbar(fig[1, 2], sc, label="log₁₀ error")
    save(joinpath(FIG, "map_legacy_error.png"), fig, px_per_unit=2)
end

function fig_history()
    lines_ = readlines(joinpath(@__DIR__, "results", "history.csv"))[2:end]
    recs = [split(l, ",") for l in lines_]
    angles = unique(r[1] for r in recs)
    impls = ["legacy+fix", "new", "new-noprune", "new-norefine"]
    for (what, col, ylabel, fname) in (("err", 4, "|c_k − c_ref|", "history_error.png"),
                                       ("points", 5, "points in spider", "history_points.png"))
        fig = Figure(size=(1200, 620))
        for (k, a) in enumerate(angles)
            ax = Axis(fig[(k-1)÷3+1, (k-1)%3+1], title="θ = $a", yscale=log10, xlabel="iteration",
                      ylabel=(k-1)%3 == 0 ? ylabel : "")
            for (j, im) in enumerate(impls)
                sel = filter(r -> r[1] == a && r[2] == im, recs)
                isempty(sel) && continue
                x = parse.(Int, getindex.(sel, 3)); y = max.(parse.(Float64, getindex.(sel, col)), 1e-18)
                lines!(ax, x, y; color=SERIES[j], linewidth=j == 1 ? 2 : max(1.5, 3.5 - 0.6j),
                       linestyle=(:solid, :solid, :dash, :dot)[j], label=im)
            end
        end
        Legend(fig[3, 1:3], fig.content[1]; orientation=:horizontal, framevisible=false)
        save(joinpath(FIG, fname), fig, px_per_unit=2)
    end
end

function fig_scaling()
    rows = read_csv(joinpath(@__DIR__, "results", "scaling.csv"))
    fams = ["cusp", "tip", "approx16", "approx14", "random", "doubling"]
    impls = ["legacy", "new-f64", "new-bigfeet", "new-bigall"]
    for (field, ylabel, fname) in ((:time, "wall time [s]", "scaling_time.png"),
                                   (:iterations, "spider iterations", "scaling_iterations.png"),
                                   (:bytes, "bytes allocated", "scaling_bytes.png"))
        fig = Figure(size=(1200, 640))
        for (k, fam) in enumerate(fams)
            ax = Axis(fig[(k-1)÷3+1, (k-1)%3+1], title=fam, yscale=log10, xlabel="period",
                      ylabel=(k-1)%3 == 0 ? ylabel : "")
            for (j, im) in enumerate(impls)
                sel = filter(r -> r.impl == im && startswith(r.note, fam * " ") && r.status != "exception", rows)
                isempty(sel) && continue
                x = [r.period for r in sel]; y = max.(Float64[getfield(r, field) for r in sel], 1e-4)
                o = sortperm(x)
                fam == "random" || lines!(ax, x[o], y[o]; color=SERIES[j], linewidth=2)
                for st in unique(r.status for r in sel)
                    s2 = [i for i in eachindex(sel) if sel[i].status == st]
                    scatter!(ax, x[s2], y[s2]; color=SERIES[j], marker=MARK[st], markersize=st == "ok" ? 7 : 11,
                             strokewidth=1, strokecolor=:white)
                end
            end
        end
        elems = [[LineElement(color=SERIES[j], linewidth=2), MarkerElement(color=SERIES[j], marker=:circle)] for j in 1:4]
        selems = [MarkerElement(color=INK, marker=MARK[s], markersize=10) for s in ("ok", "inaccurate", "unresolved", "wrong")]
        Legend(fig[3, 1:3], vcat(elems, selems), vcat(impls, ["ok", "inaccurate", "unresolved (precision)", "wrong / NaN"]);
               orientation=:horizontal, framevisible=false, nbanks=1)
        save(joinpath(FIG, fname), fig, px_per_unit=2)
    end

    # accuracy relative to component size
    fig = Figure(size=(1100, 420))
    for (k, fam) in enumerate(["tip", "approx16", "random"])
        ax = Axis(fig[1, k], title=fam, yscale=log10, xlabel="period",
                  ylabel=k == 1 ? "|c − c_ref| / component size" : "")
        for (j, im) in enumerate(impls)
            sel = filter(r -> r.impl == im && startswith(r.note, fam * " ") && isfinite(r.relerr), rows)
            isempty(sel) && continue
            scatter!(ax, [r.period for r in sel], max.([r.relerr for r in sel], 1e-80); color=SERIES[j], markersize=8, label=im)
        end
        hlines!(ax, [0.25]; color=INK, linestyle=:dash, linewidth=1)
        k == 1 && axislegend(ax; position=:lt, framevisible=false, labelsize=11)
    end
    save(joinpath(FIG, "precision.png"), fig, px_per_unit=2)
end

for f in (fig_map, fig_history, fig_scaling)
    try
        f(); println("$f ok")
    catch e
        println("$f failed: ", sprint(showerror, e))
    end
end
