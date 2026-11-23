# Run from repository root, numerical environment (no Makie required):
# julia --project=. scripts/run_constrained.jl --smoke --vacuum
# julia --project=. scripts/run_constrained.jl --standard
# julia --project=. scripts/run_constrained.jl --standard --J 30 --tmax 600 --plot
# From Julia REPL: include("scripts/run_constrained.jl"); constrained_main(["--quick"])
# Runs the constraint-solved solver, MMIEF.run_constraint_solved, and writes
# results/constrained_<mode>[_vacuum]/, including one (t, phi) file
# phi_J<J>.csv per observer.
using MMIEF
using Printf
using TOML
include(joinpath(@__DIR__, "run_options.jl"))   # --J and --tmax

# Size presets. Same domains as scripts/run_simulation.jl except --quick
# (rmax 180, tmax 75); this path uses Courant number 0.15 and KO strength 0.025.
# --vacuum sets the pulse amplitude to zero.
function constrained_config(mode::String, vacuum::Bool)
    amplitude = vacuum ? 0.0 : 1e-4
    if mode == "--smoke"
        return MMIEF.Config(rmax=70.0,dr=0.5,cfl=0.15,tmax=2.0,
               output_dt=0.25,snapshot_dt=0.5,profile_stride=2,
               pulse_amplitude=amplitude,ko_strength=0.025)
    elseif mode == "--quick"
        return MMIEF.Config(rmax=180.0,dr=0.5,cfl=0.15,tmax=75.0,
               output_dt=0.5,snapshot_dt=1.0,profile_stride=2,
               pulse_amplitude=amplitude,ko_strength=0.025)
    elseif mode == "--standard"
        return MMIEF.Config(rmax=700.0,dr=0.25,cfl=0.15,tmax=350.0,
               output_dt=0.5,snapshot_dt=2.0,profile_stride=4,
               pulse_amplitude=amplitude,ko_strength=0.025)
    elseif mode == "--long"
        return MMIEF.Config(rmax=1700.0,dr=0.25,cfl=0.15,tmax=800.0,
               output_dt=0.5,snapshot_dt=4.0,profile_stride=8,
               pulse_amplitude=amplitude,ko_strength=0.025)
    end
    error("Unknown mode $mode")
end

# The reused DataOutput module defaults to the label "partially constrained".
# Set an accurate method label in metadata.toml of a constraint-solved run.
function label_constraint_solved!(path::AbstractString)
    metapath=joinpath(path,"metadata.toml")
    meta=TOML.parsefile(metapath)
    meta["evolution"]="Hamiltonian radial mass projection and momentum A reconstruction at every RK4 stage; B evolved"
    meta["method"]="constraint_solved_mol_rk4"
    meta["core_gravity_evolved"]="B only; a constrained; A algebraic"
    open(metapath,"w") do io
        TOML.print(io,meta)
    end
    return metapath
end

# Parse flags, run the solver, save the output and relabel its metadata.
function constrained_main(args::Vector{String}=ARGS)
    if "--help" in args || "-h" in args
        println("""
        Constraint-solved MMIEF Einstein--Klein--Gordon evolution
        Options: --smoke, --quick, --standard, --long (default --quick)
                 --vacuum (zero scalar amplitude)
        Observer and run-length options:
        $RUN_OPTIONS_HELP
        Outputs: results/constrained_<mode>[_vacuum], or with --J/--tmax
                 results/constrained_<mode>[_vacuum]_T<T>_J<J>, including one
                 file phi_J<J>.csv with columns t,phi per observer
        Exactly the Hamiltonian and momentum constraints are enforced at
        each RK4 stage; B still evolves as an Einstein variable.
        See CONSTRAINED_FORMULATION.md for the terminology.
        """)
        return nothing
    end
    flags,J,tmax=parse_run_options(args)
    allowed=Set(["--smoke","--quick","--standard","--long","--vacuum","--plot"])
    all(a->a in allowed,flags) || error("Unknown option; see --help")
    modes=[a for a in flags if a in PRESET_FLAGS]
    length(modes)<=1 || error("Specify only one size")
    mode=isempty(modes) ? "--quick" : only(modes)
    vacuum="--vacuum" in flags
    p=apply_run_options(constrained_config(mode,vacuum),J,tmax)
    name=run_folder_name("constrained_"*mode[3:end]*(vacuum ? "_vacuum" : ""),p,J,tmax)
    outdir=joinpath(@__DIR__,"..","results",name)
    @printf("Running constraint-solved mode=%s, N≈%d, tmax=%g, dr=%g, observers J=%s\n",
            name,round(Int,(p.rmax-2*p.M0)/p.dr)+1,p.tmax,p.dr,string(p.observers))
    result=MMIEF.run_constraint_solved(p)
    path=MMIEF.save_run(result,outdir)
    label_constraint_solved!(path)
    @printf("Wrote %s\n",path)
    for Jobs in p.observers
        println("Saved (t, phi) at J = $Jobs: ",joinpath(path,waveform_filename(Jobs)))
    end
    "--plot" in flags && render_tail_plot(path)
    return result
end

# Run only when executed as a script, not when include()d from the REPL.
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    constrained_main(ARGS)
end
