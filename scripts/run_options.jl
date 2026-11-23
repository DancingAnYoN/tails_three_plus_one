# Command-line options shared by scripts/run_simulation.jl and
# scripts/run_constrained.jl. On top of a size preset they choose
#   --J 30            the fixed areal radius J_obs of one observer
#   --J 10,30,50      several observers (also --J=30; --J may be repeated)
#   --tmax 1000       the final coordinate time t_max (also --tmax=1000)
#   --plot            render the log–log tail plot after the run
# Values are in units of M0 (every preset uses M0 = 1). Without these options
# the preset values are used: observers J = 10, 30, 50 and the preset's tmax.
using MMIEF
using Printf

const PRESET_FLAGS = ("--smoke", "--quick", "--standard", "--long")

const RUN_OPTIONS_HELP = """
    --J J1[,J2,...]  fixed areal radii of the observers, in M0 (default 10,30,50)
    --tmax T         final coordinate time t_max, in M0 (default: preset value)
    --plot           afterwards render tail_loglog.png (log-log tail plot) with
                     scripts/render_tail.jl in the viz/ environment"""

# Parse one number given to option `name`, rejecting text, NaN and Inf.
function parse_option_number(s::AbstractString, name::AbstractString)
    x = tryparse(Float64, strip(s))
    (x !== nothing && isfinite(x)) || error("Invalid value \"$s\" for $name; expected a number")
    return x
end

"""
    parse_run_options(args) -> (flags, J, tmax)

Separate the value options `--J` and `--tmax` from the plain flags (presets,
`--vacuum`). `J` is `nothing` or the vector of requested radii; `tmax` is
`nothing` or a number. Both are in units of M0 and not yet validated.
"""
function parse_run_options(args::AbstractVector{<:AbstractString})
    flags = String[]
    J = Float64[]
    tmax = nothing
    i = 1
    while i <= length(args)
        arg = String(args[i])
        name, value = arg, nothing
        if startswith(arg, "--") && occursin('=', arg)
            name, value = split(arg, '='; limit=2)   # "--J=30" -> ("--J", "30")
        end
        if name == "--J" || name == "--tmax"
            if value === nothing                     # "--J 30": value is the next argument
                i < length(args) || error("$name needs a value; see --help")
                i += 1
                value = args[i]
            end
            if name == "--J"
                append!(J, [parse_option_number(v, name) for v in split(value, ',')])
            else
                tmax === nothing || error("Give --tmax only once")
                tmax = parse_option_number(value, name)
            end
        else
            push!(flags, arg)
        end
        i += 1
    end
    return flags, (isempty(J) ? nothing : J), tmax
end

"""
    apply_run_options(p, J, tmax) -> Config

Copy of the preset `p` with the requested observer radii `J` and final time
`tmax` (in units of M0; `nothing` keeps the preset value). Checks them against
the preset's domain before any evolution starts, then calls
`warn_outer_boundary`.
"""
function apply_run_options(p::Config, J, tmax)
    observers = J === nothing ? copy(p.observers) : p.M0 .* J
    final_t = tmax === nothing ? p.tmax : p.M0 * tmax
    final_t > 0 || error("--tmax must be positive")
    rh = 2*p.M0
    for Jobs in observers
        Jobs > rh || error(@sprintf(
            "Observer J = %g is not outside the initial horizon J_h = %g (in M0)", Jobs/p.M0, rh/p.M0))
        Jobs < p.rmax || error(@sprintf(
            "Observer J = %g is not inside this preset's outer boundary rmax = %g (in M0); use a larger preset",
            Jobs/p.M0, p.rmax/p.M0))
    end
    allunique(observers) || error("Observer radii J must be distinct")
    # Fail now rather than after the run, when save_waveforms writes the files.
    allunique(waveform_filename.(observers)) ||
        error("Observer radii J are too close: they share a phi_J<J>.csv file name")
    settings = Dict{Symbol,Any}(k => getfield(p, k) for k in fieldnames(Config))
    settings[:observers] = observers
    settings[:tmax] = final_t
    q = Config(; settings...)
    warn_outer_boundary(q)
    return q
end

"""
    run_tag(T, J) -> String

Parameter part of a run-folder name: `T600_J30`, or `T600_J10-30-50` for
several observers (values in M0, written with %g).
"""
run_tag(T::Real, J::AbstractVector{<:Real}) =
    @sprintf("T%g_J%s", T, join([@sprintf("%g", x) for x in J], "-"))

"""
    run_folder_name(base, p, J, tmax) -> String

Folder name under results/ for a preset run. Without --J and --tmax it is the
preset name `base` (e.g. "standard", "constrained_quick_vacuum"). With either
option it gets the run parameters appended (e.g. "standard_T600_J30"), so a
customised run never overwrites the preset's own results.
"""
function run_folder_name(base::AbstractString, p::Config, J, tmax)
    (J === nothing && tmax === nothing) && return String(base)
    return string(base, "_", run_tag(p.tmax/p.M0, p.observers ./ p.M0))
end

"""
    render_in_viz(script_name, run_dir)

Run scripts/<script_name> on `run_dir` in a separate Julia process with the
viz/ environment, so GLMakie is never loaded into the core environment. A
failure (e.g. viz/ not instantiated, no OpenGL) is reported, not thrown: the
data are already saved.
"""
function render_in_viz(script_name::AbstractString, run_dir::AbstractString)
    root = normpath(joinpath(@__DIR__, ".."))
    script = joinpath(root, "scripts", script_name)
    cmd = `$(Base.julia_cmd()) --project=$(joinpath(root, "viz")) $script $run_dir`
    println("Rendering with $script_name in the viz environment (GLMakie start-up takes a while)...")
    try
        run(cmd)
    catch
        # The child process has already printed its own error above.
        println("Plot failed; the data are saved. If viz/ is not set up yet, see\n",
                "README section 2, then retry with:\n",
                "  julia --project=viz scripts/$script_name ", run_dir)
    end
    return nothing
end

"""Render tail_loglog.png for `run_dir` (scripts/render_tail.jl, viz environment)."""
render_tail_plot(run_dir::AbstractString) = render_in_viz("render_tail.jl", run_dir)

"""
Warn when an observer can see the outer boundary before t_max.

In this gauge (α = a(1-β)) the ingoing scalar characteristic has coordinate
speed exactly 1. B is carried inward more slowly; on the Schwarzschild
background its speed is 1 - β(1-β) <= 1. a and A are reconstructed from the
horizon outward, and ḟ uses horizon data only. The observer's grid point
r_obs = J - f(t) moves only inward (f >= 0). Hence, for the continuum
equations, nothing from the outer boundary r = rmax reaches J before
t = rmax - J. Below that time φ(t, J) cannot contain boundary reflections. This
is a sufficient condition, not an estimate of the reflection amplitude.
"""
function warn_outer_boundary(p::Config)
    for J in p.observers
        t_causal = p.rmax - J
        if p.tmax > t_causal
            @printf("""
            WARNING: tmax/M0 = %g exceeds (rmax - J)/M0 = %g for the observer at J/M0 = %g.
                     From t/M0 = %g on, phi(t,J) is causally connected to the outer
                     boundary at rmax/M0 = %g and may contain boundary reflections.
                     Compare with a larger preset before using the late-time data.
            """, p.tmax/p.M0, t_causal/p.M0, J/p.M0, t_causal/p.M0, p.rmax/p.M0)
        end
    end
end
