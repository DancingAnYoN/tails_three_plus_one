# Run from repository root, numerical environment (no Makie required):
# julia --project=. scripts/run_constrained.jl --smoke --vacuum
# julia --project=. scripts/run_constrained.jl --standard
# From Julia REPL: include("scripts/run_constrained.jl"); constrained_main(["--quick"])
using MMIEF
using Printf
using TOML

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

function constrained_main(args::Vector{String}=ARGS)
    if "--help" in args || "-h" in args
        println("""
        Constraint-solved MMIEF Einstein--Klein--Gordon evolution
        Options: --smoke, --quick, --standard, --long (default --quick)
                 --vacuum (zero scalar amplitude)
        Outputs: results/constrained_<mode>[_vacuum]
        Exactly the Hamiltonian and momentum constraints are enforced at
        each RK4 stage; B still evolves as an Einstein variable.
        See CONSTRAINED_FORMULATION.md for the terminology.
        """)
        return nothing
    end
    allowed=Set(["--smoke","--quick","--standard","--long","--vacuum"])
    all(a->a in allowed,args) || error("Unknown option; see --help")
    modes=[a for a in args if a in ("--smoke","--quick","--standard","--long")]
    length(modes)<=1 || error("Specify only one size")
    mode=isempty(modes) ? "--quick" : only(modes)
    vacuum="--vacuum" in args
    p=constrained_config(mode,vacuum)
    name="constrained_"*mode[3:end]*(vacuum ? "_vacuum" : "")
    outdir=joinpath(@__DIR__,"..","results",name)
    @printf("Running constraint-solved mode=%s, N≈%d, tmax=%g, dr=%g\n",
            name,round(Int,(p.rmax-2*p.M0)/p.dr)+1,p.tmax,p.dr)
    result=MMIEF.run_constraint_solved(p)
    path=MMIEF.save_run(result,outdir)
    # The reused DataOutput module defaults to the label "partially
    # constrained". Set an accurate method label for these outputs.
    metapath=joinpath(path,"metadata.toml")
    meta=TOML.parsefile(metapath)
    meta["evolution"]="Hamiltonian radial mass projection and momentum A reconstruction at every RK4 stage; B evolved"
    meta["method"]="constraint_solved_mol_rk4"
    meta["core_gravity_evolved"]="B only; a constrained; A algebraic"
    open(metapath,"w") do io
        TOML.print(io,meta)
    end
    @printf("Wrote %s\n",path)
    return result
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    constrained_main(ARGS)
end
