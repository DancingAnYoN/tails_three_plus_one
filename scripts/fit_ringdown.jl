# Fundamental scalar l = 0 ringdown fit with the late-time tail included:
#     phi(T) = e^{-Γτ}(a cos ωτ + b sin ωτ) + κ tail(t),    τ = T - T1,
# where tail(t) = A t^-γ exp(B/t + C/t²) is the run's tail fit extrapolated
# back into the ringdown (tail_fit.csv from scripts/fit_tail.jl; made here from
# t >= 250 if missing), and a, b, κ, ω, Γ are fitted (MMIEF.fit_qnm_tail). For
# comparison every window is also fitted with the tail subtracted (κ = 1 fixed)
# and with no tail term. Results are compared with the Schwarzschild scalar
# l = 0 fundamental mode M_f ω = 0.110455 - 0.104896 i (exp(-iωT)), using the
# final horizon mass M_f and observer time T = t - f.
#
# Windows are given at J = 30 M0 (default [90, 160]), scanned over T1 ± 5 and
# T2 - 30, T2 + 40, and shifted to other observers by the MMIEF outgoing delay;
# windows beyond t = rmax - J (outer-boundary influence) are skipped.
#   julia --project=. scripts/fit_ringdown.jl results/T1000_J30_run [--t1 90] [--t2 160] [--no-plot]
# REPL: include("scripts/fit_ringdown.jl"); fit_ringdown("results/T1000_J30_run")
# Writes into the run directory:
#   ringdown_fit.csv        every window and method (frequencies, errors, κ, residual)
#   ringdown_fit_J<J>.csv   T, t, phi, mode, tail, model, residual for the reference window
#   ringdown_J<J>.png       via scripts/render_ringdown.jl (viz environment), unless --no-plot
# Exploratory: judge the mode by its stability over windows, observers and
# resolutions, not by one fit.
using MMIEF
using Printf
include(joinpath(@__DIR__, "fit_tail.jl"))     # fit_tail, read_tail_fits, tail_model
include(joinpath(@__DIR__, "run_options.jl"))  # render_in_viz

const QNM_WR, QNM_WI = 0.110455, 0.104896   # scalar l=0 fundamental, units of M_f
const RINGDOWN_JREF = 30.0                  # windows are specified at this J (in M0)
const RINGDOWN_METHODS = ("joint", "subtracted", "no_tail")

"""
    fit_ringdown(run_dir; t1=90.0, t2=160.0, plot=true) -> path of ringdown_fit.csv

Fit the fundamental ringdown at every observer of `run_dir` with the tail fit
included, scan the windows and write the result files (see the file header).
"""
function fit_ringdown(run_dir::AbstractString; t1::Real=90.0, t2::Real=160.0, plot::Bool=true)
    run_dir = abspath(run_dir)
    t2 > t1 + 40 || error("The reference window [t1, t2] must be longer than 40 M0")
    fits = read_tail_fits(run_dir)
    if isempty(fits)
        println("tail_fit.csv not found: fitting the tail from t = 250 first.")
        fit_tail(run_dir; tmin=250.0)
        fits = read_tail_fits(run_dir)
    end
    d = MMIEF.load_run(run_dir)
    p = d.settings; M0 = p.M0
    t = d.times; T = d.times .- d.shifts; Mf = d.masses_h[end]
    @printf("Ringdown + tail fit   run: %s   M_f = %.8f M0\n", basename(run_dir), Mf/M0)
    @printf("Reference: M_f*omega_R = %.6f, M_f*Gamma = %.6f (scalar l=0 fundamental)\n", QNM_WR, QNM_WI)
    rows = String[]
    for (j, J) in enumerate(p.observers)
        y = d.signals[j]
        key = findfirst(k -> isapprox(k, J; atol=1e-9), collect(keys(fits)))
        if key === nothing
            @printf("\nJ = %g M0: no tail fit in tail_fit.csv, skipped (run fit_tail.jl)\n", J/M0)
            continue
        end
        fp = fits[collect(keys(fits))[key]]
        template = tail_model.(Ref(fp), t)          # tail fit was made in coordinate t
        shift = mmief_outgoing_delay(J, RINGDOWN_JREF*M0, Mf)
        Tlim = T[max(searchsortedlast(t, p.rmax - J), 1)]   # boundary-free samples
        @printf("\nJ = %g M0: windows shifted by %+.2f M0 from J = %g; tail fit from t = %g (gamma = %.4f)\n",
                J/M0, shift/M0, RINGDOWN_JREF, fp.tmin, fp.gamma)
        @printf("  %-14s %11s %11s %9s %9s %9s %9s\n", "window", "Mf*omegaR", "Mf*Gamma",
                "dwR %", "dG %", "kappa", "rel rms")
        errs = Dict(m => Float64[] for m in RINGDOWN_METHODS)
        ref = nothing
        for a1 in (t1 - 5, t1, t1 + 5), a2 in (t2 - 30, t2, t2 + 40)
            # ±1e-3 M0 so that samples at the nominal ends (T = t - f, f > 0) are included.
            lo, hi = M0*a1 + shift - 1e-3*M0, M0*a2 + shift + 1e-3*M0
            if lo < T[1] || hi > Tlim
                @printf("  [%5.1f,%5.1f]  outside the boundary-free data, skipped\n", lo/M0, hi/M0)
                continue
            end
            isref = a1 == t1 && a2 == t2
            for m in RINGDOWN_METHODS
                f = m == "joint"      ? fit_qnm_tail(T, y; tlo=lo, thi=hi, tail=template, mass=Mf) :
                    m == "subtracted" ? fit_qnm_tail(T, y .- template; tlo=lo, thi=hi, mass=Mf) :
                                        fit_qnm_tail(T, y; tlo=lo, thi=hi, mass=Mf)
                f === nothing && continue
                dR = 100*(f.omegaM/QNM_WR - 1); dI = 100*(f.gammaM/QNM_WI - 1)
                append!(errs[m], (dR, dI))
                push!(rows, @sprintf("%g,%.10g,%.10g,%s,%d,%.10g,%.10g,%.6g,%.6g,%.10g,%.10g,%.10g,%.6g,%d",
                      J/M0, f.tlo/M0, f.thi/M0, m, isref, f.omegaM, f.gammaM, dR, dI,
                      f.kappa, f.a, f.b, f.rel_rms, f.n))
                if m == "joint"
                    @printf("  [%5.1f,%5.1f]%s %11.6f %11.6f %+9.3f %+9.3f %9.4f %9.2e\n",
                            f.tlo/M0, f.thi/M0, isref ? "*" : " ", f.omegaM, f.gammaM, dR, dI, f.kappa, f.rel_rms)
                    isref && (ref = f)
                end
            end
        end
        for m in RINGDOWN_METHODS
            isempty(errs[m]) && continue
            e = errs[m]
            @printf("  %-11s omegaR error %+7.2f .. %+7.2f %%   Gamma error %+7.2f .. %+7.2f %%\n", m,
                    minimum(e[1:2:end]), maximum(e[1:2:end]), minimum(e[2:2:end]), maximum(e[2:2:end]))
        end
        ref === nothing && continue
        # Decomposition of the reference-window joint fit, shown a bit beyond the window.
        w, G = ref.omegaM/Mf, ref.gammaM/Mf
        open(joinpath(run_dir, @sprintf("ringdown_fit_J%g.csv", J/M0)), "w") do io
            println(io, "T,t,phi,mode,tail,model,residual")
            for i in eachindex(T)
                (ref.tlo - 30*M0 <= T[i] <= min(ref.thi + 120*M0, Tlim)) || continue
                tau = T[i] - ref.tlo
                mode = exp(-G*tau)*(ref.a*cos(w*tau) + ref.b*sin(w*tau))
                tl = ref.kappa*template[i]
                @printf(io, "%.10g,%.10g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                        T[i], t[i], y[i], mode, tl, mode + tl, y[i] - mode - tl)
            end
        end
    end
    path = joinpath(run_dir, "ringdown_fit.csv")
    open(path, "w") do io
        println(io, "J,Tlo,Thi,method,reference,Mf_omegaR,Mf_gamma,omegaR_error_pct,gamma_error_pct,kappa,a,b,rel_rms,n")
        foreach(r -> println(io, r), rows)
    end
    println("\n* = reference window. Saved: ", path, " and ringdown_fit_J<J>.csv")
    plot && render_in_viz("render_ringdown.jl", run_dir)
    return path
end

function fit_ringdown_main(args::Vector{String}=ARGS)
    (isempty(args) || "--help" in args) &&
        return println("Usage: julia --project=. scripts/fit_ringdown.jl <run_directory> [--t1 90] [--t2 160] [--no-plot]")
    value(name, default) = (i = findfirst(==(name), args);
        i === nothing ? default : (i < length(args) ? parse(Float64, args[i+1]) : error("$name needs a value")))
    return fit_ringdown(args[1]; t1=value("--t1", 90.0), t2=value("--t2", 160.0),
                        plot=!("--no-plot" in args))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    fit_ringdown_main(ARGS)
end
