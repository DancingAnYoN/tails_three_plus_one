# Julia project runner: NO GLMakie dependency.
# Execute from project root:
#   julia --project=. scripts/run_simulation.jl --smoke --vacuum
#   julia --project=. scripts/run_simulation.jl --standard --J 30 --tmax 600 --plot
# Or REPL: include("scripts/run_simulation.jl"); main(["--smoke"])
# Runs the baseline (partially constrained) solver, MMIEF.run_evolution, and
# writes results/<mode>[_vacuum]/, including one (t, phi) file phi_J<J>.csv
# per observer. For the constraint-solved solver use scripts/run_constrained.jl.
using MMIEF
using Printf
include(joinpath(@__DIR__, "run_options.jl"))   # --J and --tmax

# Size presets. All use M0 = 1, Courant number 0.2 and KO strength 0.045;
# they differ in domain size rmax, resolution dr and run length tmax.
# --vacuum sets the pulse amplitude to zero (stationary Schwarzschild check).
function configuration(mode::String, vacuum::Bool)
    if mode=="--smoke"
        return Config(rmax=70.0,dr=0.5,cfl=0.2,tmax=2.0,output_dt=0.5,
                      snapshot_dt=0.5,profile_stride=2,ko_strength=0.045,
                      pulse_amplitude=vacuum ? 0.0 : 1e-4)
    elseif mode=="--quick"
        return Config(rmax=160.0,dr=0.5,cfl=0.2,tmax=65.0,output_dt=0.5,
                      snapshot_dt=1.0,profile_stride=2,ko_strength=0.045,
                      pulse_amplitude=vacuum ? 0.0 : 1e-4)
    elseif mode=="--standard"
        return Config(rmax=700.0,dr=0.25,cfl=0.2,tmax=350.0,output_dt=0.5,
                      snapshot_dt=2.0,profile_stride=4,ko_strength=0.045,
                      pulse_amplitude=vacuum ? 0.0 : 1e-4)
    elseif mode=="--long"
        return Config(rmax=1700.0,dr=0.25,cfl=0.2,tmax=800.0,output_dt=0.5,
                      snapshot_dt=4.0,profile_stride=8,ko_strength=0.045,
                      pulse_amplitude=vacuum ? 0.0 : 1e-4)
    else
        error("Unrecognized mode $mode")
    end
end

# Parse flags, run the solver and save the versioned output.
function main(args::Vector{String}=ARGS)
    if "--help" in args || "-h" in args
        println("""
        MMIEF Einstein–Klein–Gordon numerical simulation
        Usage: julia --project=. scripts/run_simulation.jl [options]
        --smoke | --quick | --standard | --long   (default --quick)
        --vacuum    Schwarzschild vacuum test
        $RUN_OPTIONS_HELP
        --help      This help
        Output: results/<mode>[_vacuum]/, or results/<mode>[_vacuum]_T<T>_J<J>/
        when --J or --tmax is given, containing versioned CSV + TOML and
        one file phi_J<J>.csv with columns t,phi per observer.
        Run visualization separately via --project=viz.
        """)
        return nothing
    end
    flags,J,tmax=parse_run_options(args)
    allowed=Set(["--smoke","--quick","--standard","--long","--vacuum","--plot"])
    all(x->x in allowed,flags) || error("Unrecognized flag. See --help")
    modes=[x for x in flags if x in PRESET_FLAGS]
    length(modes)<=1 || error("Specify only one simulation size")
    mode=isempty(modes) ? "--quick" : only(modes)
    p=apply_run_options(configuration(mode,"--vacuum" in flags),J,tmax)
    runname=run_folder_name(mode[3:end]*("--vacuum" in flags ? "_vacuum" : ""),p,J,tmax)
    target=joinpath(@__DIR__,"..","results",runname)
    println("Running $runname: dr=$(p.dr), tmax=$(p.tmax), rmax=$(p.rmax), observers J=$(p.observers)")
    data=run_evolution(p)
    output=save_run(data,target)
    @printf("Saved standardized output: %s\nFinal horizon mass: %.12g\n",
            output,data.masses_h[end])
    for Jobs in p.observers
        println("Saved (t, phi) at J = $Jobs: ",joinpath(output,waveform_filename(Jobs)))
    end
    "--plot" in flags && render_tail_plot(output)
    return data
end

# Run only when executed as a script, not when include()d from the REPL.
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
