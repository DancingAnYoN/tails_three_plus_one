# Fit the late-time tail at every observer of a run,
#     phi(t) = A t^(-gamma) exp(B/t + C/t^2),
# by linear least squares in ln|phi| (MMIEF.fit_tail_model), from the
# phi_J<J>.csv files (coordinate time t; for weak pulses T = t - f differs by
# < 1e-4 M0). Here B and C are just the fit coefficients of the correction
# terms, unrelated to the extrinsic curvature B = K^θ_θ and the Courant number.
#   julia --project=. scripts/fit_tail.jl results/T1000_J30_run [--tmin 150]
# REPL: include("scripts/fit_tail.jl"); fit_tail("results/T1000_J30_run")
# Prints the fit, the same fit with fewer correction terms, and its dependence
# on the window start. Writes into the run directory:
#   tail_fit.csv        one row per fit (J, window, order, A, gamma, B, C, residuals);
#                       order = number of correction terms: 0 power law, 1 + B/t, 2 + C/t^2
#   tail_fit_J<J>.csv   t, phi, phi_fit, ln_residual for the requested fit
# No GLMakie needed. The fit is exploratory: its residuals are systematic, so
# judge gamma by its window, observer and resolution dependence.
using MMIEF
using Printf
using TOML

# (J, t, phi) for every phi_J<J>.csv in run_dir, sorted by J.
function read_waveforms(run_dir::AbstractString)
    out = Tuple{Float64,Vector{Float64},Vector{Float64}}[]
    for name in readdir(run_dir)
        m = match(r"^phi_J([-+.0-9eE]+)\.csv$", name)
        m === nothing && continue
        J = parse(Float64, m.captures[1])
        lines = readlines(joinpath(run_dir, name))
        strip(lines[1]) == "t,phi" || error("Unexpected header in $name")
        rows = [parse.(Float64, split(l, ',')) for l in lines[2:end] if !isempty(strip(l))]
        push!(out, (J, [r[1] for r in rows], [r[2] for r in rows]))
    end
    isempty(out) && error("No phi_J<J>.csv files in $run_dir (runs made before they existed: rerun)")
    return sort(out; by=first)
end

const ORDER_NAMES = ("A t^-gamma", "A t^-gamma exp(B/t)", "A t^-gamma exp(B/t + C/t^2)")

"""
    read_tail_fits(run_dir) -> Dict(J => parameters)

Fits written by `fit_tail` to tail_fit.csv: per observer, the first full-model
row (order 2), i.e. the fit on the requested window, as
(tmin, tmax, A, gamma, B, C, rms). Empty if the file does not exist.
"""
function read_tail_fits(run_dir::AbstractString)
    fits = Dict{Float64,NamedTuple}()
    path = joinpath(run_dir, "tail_fit.csv")
    isfile(path) || return fits
    lines = readlines(path)
    header = String.(strip.(split(lines[1], ',')))
    col(name) = something(findfirst(==(name), header), 0)
    for line in lines[2:end]
        isempty(strip(line)) && continue
        v = split(line, ',')
        num(name) = parse(Float64, v[col(name)])
        J = num("J")
        (num("order") == 2 && !haskey(fits, J)) || continue
        fits[J] = (tmin=num("tmin"), tmax=num("tmax"), A=num("A"), gamma=num("gamma"),
                   B=num("B"), C=num("C"), rms=num("rms_ln_residual"))
    end
    return fits
end

"""Value of a fitted tail A t^-γ exp(B/t + C/t²) (parameters from `read_tail_fits`) at t."""
tail_model(fp, t::Real) = fp.A * t^(-fp.gamma) * exp(fp.B/t + fp.C/t^2)

"""
    fit_tail(run_dir; tmin=150.0) -> path of tail_fit.csv

Fit every observer of `run_dir` on tmin <= t <= min(t_end, rmax - J). The upper
limit keeps the window free of outer-boundary influence.
"""
function fit_tail(run_dir::AbstractString; tmin::Real=150.0)
    run_dir = abspath(run_dir)
    cfg = TOML.parsefile(joinpath(run_dir, "config.toml"))
    M0 = Float64(cfg["M0"]); rmax = Float64(cfg["rmax"])
    rows = String[]
    println("Tail fit  phi(t) = A t^-gamma exp(B/t + C/t^2)   run: ", basename(run_dir))
    println("Times in M0 (M0 = $M0); residuals are in ln|phi| (relative residuals of phi).")
    for (J, t, phi) in read_waveforms(run_dir)
        tmax = min(t[end], rmax - J)    # boundary-free part (scripts/run_options.jl)
        tmax > tmin + 50 || (println("\nJ = $J: window [$tmin, $tmax] too short, skipped"); continue)
        @printf("\nJ = %g M0, window t in [%g, %g]%s\n", J, tmin, tmax,
                tmax < t[end] ? @sprintf(" (cut at rmax - J; data end at %g)", t[end]) : "")
        @printf("  %-30s %12s %9s %10s %11s %10s %10s\n",
                "model", "A", "gamma", "B", "C", "rms res", "max res")
        fits = try
            [fit_tail_model(t, phi; tmin=tmin, tmax=tmax, order=k) for k in 0:2]
        catch err
            err isa ErrorException || rethrow()
            println("  skipped: ", err.msg)    # e.g. phi still changes sign: use a later --tmin
            continue
        end
        for (k, f) in enumerate(fits)         # k = order + 1
            @printf("  %-30s %12.5e %9.5f %10.3f %11.2f %10.2e %10.2e\n",
                    ORDER_NAMES[k], f.A, f.gamma, f.b1, f.b2, f.rms, f.maxabs)
            push!(rows, @sprintf("%g,%g,%g,%d,%.12g,%.12g,%.12g,%.12g,%.6g,%.6g,%d",
                  J, f.tmin, f.tmax, k - 1, f.A, f.gamma, f.b1, f.b2, f.rms, f.maxabs, f.n))
        end
        # Window-start dependence of the full model (systematic uncertainty).
        println("  window-start dependence of the full model:")
        @printf("  %8s %9s %10s %11s %10s\n", "t_min", "gamma", "B", "C", "rms res")
        for t1 in unique(sort([Float64(tmin), 200.0, 250.0, 300.0, 400.0, 500.0]))
            (t1 >= tmin && t1 <= tmax - 150) || continue
            f = fit_tail_model(t, phi; tmin=t1, tmax=tmax, order=2)
            @printf("  %8g %9.5f %10.3f %11.2f %10.2e\n", t1, f.gamma, f.b1, f.b2, f.rms)
            t1 == tmin || push!(rows, @sprintf("%g,%g,%g,%d,%.12g,%.12g,%.12g,%.12g,%.6g,%.6g,%d",
                  J, f.tmin, f.tmax, 2, f.A, f.gamma, f.b1, f.b2, f.rms, f.maxabs, f.n))
        end
        # Data, fit and residual of the requested fit.
        f = fits[3]
        open(joinpath(run_dir, @sprintf("tail_fit_J%g.csv", J)), "w") do io
            println(io, "t,phi,phi_fit,ln_residual")
            for i in eachindex(t)
                tmin <= t[i] <= tmax || continue
                model = f.A * t[i]^(-f.gamma) * exp(f.b1/t[i] + f.b2/t[i]^2)
                @printf(io, "%.17g,%.17g,%.17g,%.6e\n", t[i], phi[i], model, log(phi[i]/model))
            end
        end
    end
    path = joinpath(run_dir, "tail_fit.csv")
    open(path, "w") do io
        println(io, "J,tmin,tmax,order,A,gamma,B,C,rms_ln_residual,max_ln_residual,n")
        foreach(r -> println(io, r), rows)
    end
    println("\nSaved: ", path, " and tail_fit_J<J>.csv (data, fit, residual)")
    return path
end

function fit_tail_main(args::Vector{String}=ARGS)
    (isempty(args) || "--help" in args) &&
        return println("Usage: julia --project=. scripts/fit_tail.jl <run_directory> [--tmin 150]")
    tmin = 150.0
    i = findfirst(==("--tmin"), args)
    if i !== nothing
        i < length(args) || error("--tmin needs a value")
        tmin = parse(Float64, args[i+1])
    end
    return fit_tail(args[1]; tmin=tmin)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    fit_tail_main(ARGS)
end
