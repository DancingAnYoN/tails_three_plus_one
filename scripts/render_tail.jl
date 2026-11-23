# Late-time tail plot for a saved run: |φ| against observer time T = t - f on
# log–log axes at every fixed-areal-radius observer, with a T⁻³ reference line,
# and the local power index p(T) = -d ln|φ| / d ln T.
#   julia --project=viz scripts/render_tail.jl results/standard [--show]
# Writes tail_loglog.png into the run directory and prints a short p(T)
# summary. The run scripts call this with --plot. Reads the saved output only;
# the solver is not rerun.
using MMIEF
using GLMakie
using Printf
include(joinpath(@__DIR__, "fit_tail.jl"))   # read_tail_fits, tail_model

const LPI_HALFWIDTH = 25.0   # half-width of the local log–log fit window, in M0
const PALETTE = [:dodgerblue, :darkorange, :seagreen, :purple, :firebrick]

"""
    make_tail_plot(run_dir; show_window=false) -> path

Save `tail_loglog.png` in `run_dir`. Left: log10|φ| against log10 T. Right:
p(T) from `local_loglog_slope`, defined only where φ keeps one sign over the fit
window. Dotted lines mark t = rmax - J, after which the outer boundary can
influence the observer.
"""
function make_tail_plot(run_dir::AbstractString; show_window::Bool=false)
    data = MMIEF.load_run(run_dir)
    p = data.settings
    M = p.M0
    T = data.times .- data.shifts                    # observer time T = t - f
    if !any(s -> any(!iszero, s), data.signals)
        println("phi is identically zero at every observer (vacuum run?): no tail to plot.")
        return nothing
    end
    GLMakie.activate!()

    fig = Figure(size=(1500, 720), fontsize=17)
    Label(fig[0, 1:2], @sprintf("Late-time scalar tail | %s | tmax = %g M₀, rmax = %g M₀",
                                basename(normpath(run_dir)), p.tmax/M, p.rmax/M),
          fontsize=22, font=:bold)
    ax1 = Axis(fig[1, 1], xscale=log10, yscale=log10,
               xlabel="T / M₀", ylabel="|φ(T, J)|", title="Log–log amplitude")
    ax2 = Axis(fig[1, 2], xlabel="T / M₀", ylabel="p(T) = −d ln|φ| / d ln T",
               title=@sprintf("Local power index (fit window ±%g M₀)", LPI_HALFWIDTH))

    # y-range: from below the late-time amplitudes to above the peak; the deep
    # spikes at zero crossings of the ringdown are clipped.
    ymax = maximum(maximum(abs, s) for s in data.signals)
    ylate = minimum(minimum(abs, s[max(1, end - length(s) ÷ 10):end]) for s in data.signals)
    ylo = max(0.2*ylate, 1e-16*ymax)

    fits = read_tail_fits(run_dir)
    println("Local power index p(T) = -d ln|phi| / d ln T (exploratory):")
    guide = nothing        # (T range, amplitude at its end) for the T⁻³ line
    any_index = false      # some p(T) value lies inside the plotted range
    for (k, J) in enumerate(p.observers)
        y = data.signals[k]
        col = PALETTE[mod1(k, length(PALETTE))]
        name = @sprintf("J = %g M₀", J/M)
        keep = findall(i -> T[i] > 0 && y[i] != 0, eachindex(T))
        lines!(ax1, T[keep] ./ M, abs.(y[keep]); color=col, linewidth=2, label=name)

        # Last sample not yet causally connected to the outer boundary:
        # coordinate time t <= rmax - J (see scripts/run_options.jl).
        iclean = searchsortedlast(data.times, p.rmax - J)
        if iclean < length(T)
            for ax in (ax1, ax2)
                vlines!(ax, [T[max(iclean, 1)]/M]; color=col, linestyle=:dot, linewidth=2)
            end
        end
        iclean = max(iclean, 1)
        # Tail region: within the boundary-free samples, after both the peak of
        # |φ| and the last sign change of φ (a sign change after iclean may be
        # a reflection and does not move it). Earlier sign-definite stretches
        # (rising pulse, between ringdown zero crossings) do not measure a tail.
        zc = findlast(i -> sign(y[i]) != sign(y[i+1]), 1:iclean-1)
        ipeak = argmax(abs.(view(y, 1:iclean)))
        itail = max(zc === nothing ? 1 : zc + 1, ipeak)

        # p(T) only in the tail region, with fit windows cut off at its start.
        slope = local_loglog_slope(T, y; halfwidth=LPI_HALFWIDTH*M, tstart=T[itail])
        lines!(ax2, T ./ M, slope; color=col, linewidth=2, label=name)
        any_index |= any(x -> isfinite(x) && 0 <= x <= 8, slope)

        # T⁻³ guide for the first observer: through its last boundary-free
        # sample, over the tail region (at most a decade).
        if k == 1 && y[iclean] != 0 && T[iclean] > 0
            Ta = T[iclean]
            Tb = max(Ta/10, T[itail])
            Tb < Ta && (guide = (range(Tb, Ta, length=50), abs(y[iclean])))
        end

        # Fitted model from scripts/fit_tail.jl in red on its window: the curve,
        # and its local exponent p_fit = γ + B/t + 2C/t² (fit in coordinate t;
        # T = t - f differs by < 1e-4 M0 for weak pulses).
        kf = findfirst(x -> isapprox(x, J; atol=1e-9), collect(keys(fits)))
        if kf !== nothing
            fp = fits[collect(keys(fits))[kf]]
            tf = range(fp.tmin, fp.tmax, length=400)
            flabel = @sprintf("fit, J = %g M₀, t ≥ %g M₀: γ = %.4f", J/M, fp.tmin/M, fp.gamma)
            lines!(ax1, tf ./ M, abs.(tail_model.(Ref(fp), tf));
                   color=:red, linestyle=:dash, linewidth=3, label=flabel)
            lines!(ax2, tf ./ M, fp.gamma .+ fp.B ./ tf .+ 2 .* fp.C ./ tf.^2;
                   color=:red, linestyle=:dash, linewidth=3, label=flabel)
            @printf("  J = %g M0: fit A t^-gamma exp(B/t + C/t^2) on t in [%g, %g]: A = %.5e, gamma = %.5f, B = %.3f, C = %.2f, rms ln-residual = %.2e\n",
                    J/M, fp.tmin, fp.tmax, fp.A, fp.gamma, fp.B, fp.C, fp.rms)
        end

        # Console summary over the last quarter of the boundary-free tail region.
        good = [i for i in itail:iclean if isfinite(slope[i])]
        if isempty(good)
            @printf("  J = %g M0: no tail window yet (peak of |phi| at T = %.1f M0, last sign change at T = %s M0)\n",
                    J/M, T[ipeak]/M, zc === nothing ? "none" : @sprintf("%.1f", T[zc]/M))
        else
            Tlo = T[good[1]] + 0.75*(T[good[end]] - T[good[1]])
            last_q = [slope[i] for i in good if T[i] >= Tlo]
            @printf("  J = %g M0: candidate tail T >= %.1f M0; p(T = %.1f M0) = %.4f; p over T in [%.1f, %.1f] M0: %.4f .. %.4f\n",
                    J/M, T[itail]/M, T[good[end]]/M, slope[good[end]],
                    Tlo/M, T[good[end]]/M, minimum(last_q), maximum(last_q))
        end
    end
    if guide !== nothing
        Tg, ya = guide
        lines!(ax1, Tg ./ M, ya .* (Tg[end] ./ Tg).^3;
               color=:black, linestyle=:dash, linewidth=2, label="∝ T⁻³ (l = 0 Price tail)")
    end
    hlines!(ax2, [3.0]; color=:black, linestyle=:dash, linewidth=2, label="p = 3")
    any_index || text!(ax2, 0.5, 0.55; space=:relative, align=(:center, :center),
                       text="No sign-definite tail window in range yet:\nincrease --tmax",
                       fontsize=16, color=:gray40)

    xlims!(ax1, max(1.0, T[findfirst(>(0), T)]/M), 1.1*T[end]/M)
    ylims!(ax1, ylo, 3*ymax)
    xlims!(ax2, 0, T[end]/M)
    ylims!(ax2, 0, 8)
    axislegend(ax1, position=:lb)
    axislegend(ax2, position=:rt)
    Label(fig[2, 1:2],
          "T = t − f is observer time. p(T) is a local least-squares slope of ln|φ| against ln T, shown after the peak and the last sign change of φ, where φ keeps one sign over the window.\n" *
          "Dotted lines: t = rmax − J; later samples can contain outer-boundary reflections. Exploratory: compare resolutions, observers and rmax before quoting p.",
          fontsize=13, tellwidth=false, justification=:left)

    outpath = joinpath(abspath(run_dir), "tail_loglog.png")
    save(outpath, fig)
    println("Saved tail plot: ", outpath)
    show_window && display(fig)
    return outpath
end

function tail_main(args::Vector{String}=ARGS)
    isempty(args) && error("Usage: render_tail.jl <run_directory> [--show]")
    all(x -> x == "--show", args[2:end]) || error("Unknown render option")
    return make_tail_plot(args[1]; show_window=("--show" in args))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    tail_main(ARGS)
end
