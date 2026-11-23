# One-command constraint-solved run for a total time T and observer radius J
# (both in units of M0). From the project root:
#   julia --project=. scripts/run_TJ.jl 600 30
#   julia --project=. scripts/run_TJ.jl 600 10,30,50       several observers
#   julia --project=. scripts/run_TJ.jl 600 30 --dr 0.5    coarser grid
#   julia --project=. scripts/run_TJ.jl 600 30 --no-plot   skip the tail plot
# REPL: include("scripts/run_TJ.jl"); run_TJ(600, 30)
#
# Solver: MMIEF.run_constraint_solved. a and A come from the Hamiltonian and
# momentum constraints at every RK4 stage, and B is evolved, so this is not
# strictly a fully constrained scheme (see CONSTRAINED_FORMULATION.md).
# The outer boundary is placed at rmax >= T + J + 10 M0, so no outer-boundary
# signal can reach the observer during the run (see scripts/run_options.jl).
# Output: results/T<T>_J<J>_run/, e.g. results/T600_J30_run/, holding the
# standard run files, phi_J<J>.csv with (t, phi) and tail_loglog.png.
include(joinpath(@__DIR__, "run_constrained.jl"))   # also loads run_options.jl

const TJ_DR = 0.25        # default radial spacing in M0 (as the --standard preset)
const TJ_CFL = 0.125      # dt = dr/8: for dr = 2^-k the samples fall exactly on multiples of output_dt
const TJ_MARGIN = 10.0    # outer boundary at least this far beyond T + J, in M0
const TJ_MIN_RMAX = 70.0  # keeps the initial pulse (r0 = 20, σ = 2) far from the boundary

"""
    tj_run_name(T, J; dr=TJ_DR) -> String

Folder name of a run: `T600_J30_run`; several observers `T600_J10-30-50_run`;
a non-default grid adds the spacing, `T600_J30_dr0.5_run`.
"""
function tj_run_name(T::Real, J::AbstractVector{<:Real}; dr::Real=TJ_DR)
    name = run_tag(T, J)
    dr == TJ_DR || (name *= @sprintf("_dr%g", dr))
    return name * "_run"
end

"""
    tj_config(T, J; dr=TJ_DR) -> Config

Constraint-solved settings (pulse, KO strength 0.025, output_dt 0.5 from the
--standard preset) with final time T, observers J and the smallest
grid-aligned outer boundary rmax >= max(T + max(J) + TJ_MARGIN, TJ_MIN_RMAX).
Snapshot cadence and stride grow with the run, so snapshots.csv holds at most
about 300 frames of about 1000 points.
"""
function tj_config(T::Real, J::AbstractVector{<:Real}; dr::Real=TJ_DR)
    (isfinite(T) && T > 0) || error("T must be a positive number")
    (isfinite(dr) && dr > 0) || error("--dr must be a positive number")
    isempty(J) && error("Give at least one observer radius J")
    base = constrained_config("--standard", false)
    rh = 2*base.M0
    rmin = base.M0*max(T + maximum(J) + TJ_MARGIN, TJ_MIN_RMAX)
    rmax = rh + dr*ceil((rmin - rh)/dr - 1e-9)      # (rmax - 2M0)/dr is an integer
    N = round(Int, (rmax - rh)/dr) + 1
    settings = Dict{Symbol,Any}(k => getfield(base, k) for k in fieldnames(Config))
    settings[:rmax] = rmax
    settings[:dr] = Float64(dr)
    settings[:cfl] = TJ_CFL
    settings[:snapshot_dt] = max(2.0, ceil(T/300))
    settings[:profile_stride] = max(4, cld(N, 1000))
    # Validates J and sets observers and tmax; no warning can trigger here.
    return apply_run_options(Config(; settings...), collect(Float64, J), Float64(T))
end

"""
    run_TJ(T, J; dr=0.25, plot=true) -> RunData

Run the constraint-solved scheme to time T with observer(s) at areal radius J
(a number or a vector, in M0), save to results/<tj_run_name>/ and, if `plot`,
render tail_loglog.png there. If the folder exists, the run files (and
tail_loglog.png when plotting) are overwritten; other files in it are kept.
"""
function run_TJ(T::Real, J; dr::Real=TJ_DR, plot::Bool=true)
    Js = J isa Real ? [Float64(J)] : collect(Float64, J)
    p = tj_config(T, Js; dr=dr)
    outdir = normpath(joinpath(@__DIR__, "..", "results", tj_run_name(T, Js; dr=dr)))
    N = round(Int, (p.rmax - 2*p.M0)/p.dr) + 1
    steps = ceil(Int, p.tmax/(p.cfl*p.dr))
    @printf("Constraint-solved run: T = %g M0, J = %s M0\n", T, join([@sprintf("%g", x) for x in Js], ", "))
    @printf("  grid: rmax = %g M0 (>= T + J + %g), dr = %g, N = %d; %d time steps of %g M0\n",
            p.rmax, TJ_MARGIN, p.dr, N, steps, p.cfl*p.dr)
    @printf("  estimated solver time: %.0f s (about 2 µs per grid point and step)\n", 2e-6*N*steps)
    isdir(outdir) && println("  overwriting the run files in existing ", outdir,
                             " (other files there are kept)")
    result = MMIEF.run_constraint_solved(p)
    path = MMIEF.save_run(result, outdir)
    label_constraint_solved!(path)
    println("Saved run: ", path)
    for Jobs in p.observers
        println("Saved (t, phi) at J = $Jobs: ", joinpath(path, waveform_filename(Jobs)))
    end
    plot && render_tail_plot(path)
    return result
end

const TJ_HELP = """
Constraint-solved MMIEF run for total time T and areal radius J (units of M0).
Usage: julia --project=. scripts/run_TJ.jl T J [--dr DR] [--no-plot]
  T          final time, e.g. 600
  J          observer areal radius, e.g. 30, or a list 10,30,50
  --dr DR    radial spacing (default 0.25); also appended to the folder name
  --no-plot  do not render tail_loglog.png (it needs the viz/ environment)
Output: results/T<T>_J<J>_run/, e.g. results/T600_J30_run/phi_J30.csv.
The outer boundary is placed at rmax >= T + J + 10, beyond causal reach."""

# Command-line front end: two positional numbers T and J plus options.
function run_TJ_main(args::Vector{String}=ARGS)
    if isempty(args) || "--help" in args || "-h" in args
        println(TJ_HELP)
        return nothing
    end
    positional = String[]
    dr = TJ_DR
    plot = true
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--no-plot"
            plot = false
        elseif a == "--dr"
            i < length(args) || error("--dr needs a value; see --help")
            i += 1
            dr = parse_option_number(args[i], "--dr")
        elseif startswith(a, "--dr=")
            dr = parse_option_number(a[6:end], "--dr")
        elseif startswith(a, "--")
            error("Unknown option $a; see --help")
        else
            push!(positional, a)
        end
        i += 1
    end
    length(positional) == 2 || error("Give exactly two values, T and J; see --help")
    T = parse_option_number(positional[1], "T")
    J = [parse_option_number(v, "J") for v in split(positional[2], ',')]
    return run_TJ(T, J; dr=dr, plot=plot)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_TJ_main(ARGS)
end
