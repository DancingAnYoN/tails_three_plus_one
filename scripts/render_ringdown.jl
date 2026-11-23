# Visualise the ringdown + tail fit of scripts/fit_ringdown.jl: one figure
# ringdown_J<J>.png per observer, from ringdown_fit.csv and ringdown_fit_J<J>.csv.
#   julia --project=viz scripts/render_ringdown.jl results/T1000_J30_run [--show]
# Panels: signal with the joint model and its two parts; the ringdown with the
# fitted tail removed; the residual; the frequencies of all windows and methods
# against the Schwarzschild scalar l = 0 fundamental mode.
using GLMakie
using Printf

const RD_WR, RD_WI = 0.110455, 0.104896   # reference M_f ω_R, M_f Γ

# CSV with a header row -> Dict(column => vector of strings).
function read_columns(path::AbstractString)
    lines = filter(l -> !isempty(strip(l)), readlines(path))
    header = String.(strip.(split(lines[1], ',')))
    cols = Dict(h => String[] for h in header)
    for line in lines[2:end], (h, v) in zip(header, split(line, ','))
        push!(cols[h], String(strip(v)))
    end
    return cols
end
num(cols, name) = parse.(Float64, cols[name])

"""
    make_ringdown_plot(run_dir, J; show_window=false) -> path

Save ringdown_J<J>.png in `run_dir` for the observer at areal radius J (in M0).
"""
function make_ringdown_plot(run_dir::AbstractString, J::Real; show_window::Bool=false)
    dec = read_columns(joinpath(run_dir, @sprintf("ringdown_fit_J%g.csv", J)))
    T = num(dec, "T"); phi = num(dec, "phi"); mode = num(dec, "mode")
    tl = num(dec, "tail"); model = num(dec, "model"); res = num(dec, "residual")
    all_fits = read_columns(joinpath(run_dir, "ringdown_fit.csv"))
    sel = findall(i -> isapprox(parse(Float64, all_fits["J"][i]), J; atol=1e-9),
                  eachindex(all_fits["J"]))
    method = all_fits["method"][sel]; isref = all_fits["reference"][sel] .== "1"
    wR = num(all_fits, "Mf_omegaR")[sel]; G = num(all_fits, "Mf_gamma")[sel]
    Tlo = num(all_fits, "Tlo")[sel]; Thi = num(all_fits, "Thi")[sel]
    kap = num(all_fits, "kappa")[sel]
    r = findfirst(i -> isref[i] && method[i] == "joint", eachindex(method))
    r === nothing && error("No reference-window joint fit for J = $J in ringdown_fit.csv")

    GLMakie.activate!()
    fig = Figure(size=(1500, 1060), fontsize=17)
    Label(fig[0, 1:2], @sprintf("Ringdown + tail fit | %s | J = %g M₀ | window [%.0f, %.0f] M₀: M_f ω = %.5f − %.5f i (κ = %.3f)",
                                basename(normpath(run_dir)), J, Tlo[r], Thi[r], wR[r], G[r], kap[r]),
          fontsize=21, font=:bold)
    pos(x) = [abs(v) > 0 ? abs(v) : NaN for v in x]   # |x| for log axes, zeros hidden
    shade!(ax) = vspan!(ax, [Tlo[r]], [Thi[r]]; color=(:gray, 0.15))

    ax1 = Axis(fig[1, 1], yscale=log10, xlabel="T / M₀", ylabel="|φ|",
               title="Signal and joint model  φ = mode + κ·tail")
    shade!(ax1)
    lines!(ax1, T, pos(phi); color=:dodgerblue, linewidth=3, label="data")
    lines!(ax1, T, pos(model); color=:red, linestyle=:dash, linewidth=2.5, label="mode + κ·tail")
    lines!(ax1, T, pos(tl); color=:black, linestyle=:dot, linewidth=2, label="κ·tail (extrapolated fit)")
    lines!(ax1, T, pos(mode); color=:darkorange, linestyle=:dashdot, linewidth=2, label="fitted mode")
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[1, 2], yscale=log10, xlabel="T / M₀", ylabel="|φ − κ·tail|",
               title="Ringdown with the fitted tail removed")
    shade!(ax2)
    lines!(ax2, T, pos(phi .- tl); color=:dodgerblue, linewidth=3, label="data − κ·tail")
    lines!(ax2, T, pos(mode); color=:red, linestyle=:dash, linewidth=2.5,
           label=@sprintf("mode, M_f ω_R = %.4f, M_f Γ = %.4f", wR[r], G[r]))
    # Envelope sqrt(a² + b²) e^{-Γτ} of the fitted mode.
    a = num(all_fits, "a")[sel][r]; b = num(all_fits, "b")[sel][r]
    lines!(ax2, T, hypot(a, b) .* exp.(-G[r] .* (T .- Tlo[r])); color=:black, linestyle=:dot,
           linewidth=2, label="mode envelope ∝ exp(−Γτ)")
    axislegend(ax2, position=:rt)

    ax3 = Axis(fig[2, 1], yscale=log10, xlabel="T / M₀", ylabel="|residual|",
               title="Residual of the joint fit  |φ − model|")
    shade!(ax3)
    lines!(ax3, T, pos(phi); color=(:gray, 0.6), linewidth=2, label="|φ| for scale")
    lines!(ax3, T, pos(res); color=:dodgerblue, linewidth=2, label="|φ − model|")
    axislegend(ax3, position=:rt)

    # Frequencies of all windows; zoomed to ±8 % of the reference, off-scale fits counted.
    ax4 = Axis(fig[2, 2], xlabel="M_f ω_R", ylabel="M_f Γ",
               title="All windows and methods vs Schwarzschild l = 0 fundamental")
    style = Dict("joint" => (:red, :circle, "joint: mode + κ·tail"),
                 "subtracted" => (:dodgerblue, :rect, "tail subtracted (κ = 1)"),
                 "no_tail" => (:gray50, :utriangle, "no tail term"))
    xl = (0.92RD_WR, 1.08RD_WR); yl = (0.92RD_WI, 1.08RD_WI)
    offscale = String[]
    for m in ("no_tail", "subtracted", "joint")
        k = findall(==(m), method)
        isempty(k) && continue
        inside = [i for i in k if xl[1] <= wR[i] <= xl[2] && yl[1] <= G[i] <= yl[2]]
        length(inside) < length(k) && push!(offscale, @sprintf("%s: %d of %d off scale", m, length(k) - length(inside), length(k)))
        col, mk, lab = style[m]
        scatter!(ax4, wR[inside], G[inside]; color=col, marker=mk, markersize=14, label=lab)
    end
    scatter!(ax4, [wR[r]], [G[r]]; color=:transparent, strokecolor=:red, strokewidth=3,
             markersize=26, label="reference window (joint)")
    scatter!(ax4, [RD_WR], [RD_WI]; color=:black, marker=:star5, markersize=24,
             label="Schwarzschild 0.110455 − 0.104896 i")
    xlims!(ax4, xl...); ylims!(ax4, yl...)
    axislegend(ax4, position=:lb)
    isempty(offscale) || text!(ax4, 0.98, 0.97; space=:relative, align=(:right, :top),
                               text=join(offscale, "\n"), fontsize=14, color=:gray30)

    Label(fig[3, 1:2],
          "Shaded: reference fit window. The tail term is the t ≥ tmin tail fit A t^−γ exp(B/t + C/t²) extrapolated back, scaled by a fitted κ. " *
          "Frequencies in units of the final horizon mass M_f; T = t − f.\n" *
          "Exploratory single-mode fits: compare windows, observers and resolutions before identifying a QNM.",
          fontsize=13, tellwidth=false, justification=:left)

    path = joinpath(abspath(run_dir), @sprintf("ringdown_J%g.png", J))
    save(path, fig)
    println("Saved ringdown plot: ", path)
    show_window && display(fig)
    return path
end

# Plot every observer that has a ringdown_fit_J<J>.csv in run_dir.
function render_ringdown_main(args::Vector{String}=ARGS)
    isempty(args) && error("Usage: render_ringdown.jl <run_directory> [--show]")
    all(x -> x == "--show", args[2:end]) || error("Unknown render option")
    run_dir = args[1]
    Js = Float64[]
    for name in readdir(run_dir)
        m = match(r"^ringdown_fit_J([-+.0-9eE]+)\.csv$", name)
        m === nothing || push!(Js, parse(Float64, m.captures[1]))
    end
    isempty(Js) && error("No ringdown_fit_J<J>.csv in $run_dir: run scripts/fit_ringdown.jl first")
    return [make_ringdown_plot(run_dir, J; show_window=("--show" in args)) for J in sort(Js)]
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    render_ringdown_main(ARGS)
end
