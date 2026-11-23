# In-memory record of a run (RunData) and its portable on-disk form:
# config.toml, metadata.toml, diagnostics.csv, observers.csv, snapshots.csv
# and final_state.csv (schema version SCHEMA_VERSION), plus one two-column
# waveform file phi_J<J>.csv per observer. Uses only Julia stdlib, so saved
# runs can be loaded by the separate visualization environment.
module DataOutput
using Printf, TOML
using ..Types: Config, IA, IB, IPHI, IPI, IFIELD
using ..Diagnostics: diagnostics
using ..Observers: sample_fixed_J
export RunData, store_sample!, store_frame!, save_csv, save_run, load_run,
       save_waveforms, waveform_filename

# Bump when the file layout changes; load_run rejects other versions.
const SCHEMA_VERSION = 1
mutable struct RunData
    settings::Config
    r::Vector{Float64}                   # fixed grid coordinates
    # Time series, one entry per sample (every output_dt):
    times::Vector{Float64}               # coordinate time t
    signals::Vector{Vector{Float64}}     # φ at each fixed J_obs (one vector per observer)
    shifts::Vector{Float64}              # f(t); observer time is T = t - f
    masses_h::Vector{Float64}            # horizon mass J_h/2
    masses_out::Vector{Float64}          # Misner–Sharp mass at r_max
    ham_l2::Vector{Float64}              # RMS Hamiltonian residual
    scalar_constraint_l2::Vector{Float64} # RMS of Φ - ∂_r φ
    horizon_error::Vector{Float64}       # aJB - 1 at r_h after projection
    horizon_preprojection::Vector{Float64} # aJB - 1 at r_h before projection
    # Radial snapshots, one entry per frame (every snapshot_dt):
    frame_times::Vector{Float64}
    frame_shift::Vector{Float64}
    frame_profiles::Vector{Vector{Float32}} # Jφ at every profile_stride-th grid point
    final_state::Matrix{Float64}         # state u at the end (no t or f: not a restart file)
end

"""Empty `RunData` for configuration `p`, grid `r` and initial state `u`."""
function RunData(p::Config,r::Vector{Float64},u::Matrix{Float64})
    return RunData(p,r,Float64[],[Float64[] for _ in p.observers],Float64[],
                   Float64[],Float64[],Float64[],Float64[],Float64[],
                   Float64[],Float64[],Float64[],Vector{Float32}[],u)
end

"""
Append one time-series sample: φ at every observer plus the diagnostics of
`(u, f)` at time `t`. `preprojection_error` is the horizon residual returned by
the RK4 step before projection.
"""
function store_sample!(data::RunData,t::Float64,u::Matrix{Float64},f::Float64,
                       preprojection_error::Float64=0.0)
    # Sample all observers before modifying RunData, so a failed interpolation
    # cannot leave its arrays with inconsistent lengths.
    wave = [sample_fixed_J(u,data.r,f,J,data.settings;time=t)
            for J in data.settings.observers]
    d=diagnostics(u,f,data.r,data.settings)
    push!(data.times,t)
    push!(data.shifts,f)
    push!(data.masses_h,d.M_h)
    push!(data.masses_out,d.M_out)
    push!(data.ham_l2,d.ham_l2)
    push!(data.scalar_constraint_l2,d.scalar_constraint_l2)
    push!(data.horizon_error,d.horizon_error)
    push!(data.horizon_preprojection,preprojection_error)
    for k in eachindex(data.settings.observers)
        push!(data.signals[k],wave[k])
    end
end

"""
Append one radial snapshot of Jφ (the rescaled field, roughly constant along
outgoing rays far away) on every `profile_stride`-th grid point, in Float32.
"""
function store_frame!(data::RunData,t::Float64,u::Matrix{Float64},f::Float64)
    push!(data.frame_times,t)
    push!(data.frame_shift,f)
    stride=data.settings.profile_stride
    push!(data.frame_profiles, Float32[
        (data.r[i]+f)*u[i,IFIELD] for i in 1:stride:length(data.r)])
end
"""
Write all time-series samples to one wide CSV file (one column per observer).
Legacy single-file format; `save_run` is the versioned format.
"""
function save_csv(data::RunData,filename::AbstractString)
    open(filename,"w") do io
        cols=["t","T_asymptotic","f","M_h","M_outer","hamiltonian_L2","scalar_definition_L2","horizon_error",
              "horizon_preprojection_error"]
        for J in data.settings.observers
            push!(cols,@sprintf("phi_J%.2f",J))
        end
        println(io,join(cols,","))
        for i in eachindex(data.times)
            vals=[data.times[i],data.times[i]-data.shifts[i],data.shifts[i],data.masses_h[i],data.masses_out[i],
                  data.ham_l2[i],data.scalar_constraint_l2[i],
                  data.horizon_error[i],data.horizon_preprojection[i]]
            append!(vals,[sig[i] for sig in data.signals])
            println(io,join([@sprintf("%.16e",x) for x in vals],","))
        end
    end
    return filename
end

"""
    waveform_filename(J_obs) -> String

Name of the two-column waveform file for observer radius `J_obs`, written with
`%g`: `phi_J30.csv` for J_obs = 30, `phi_J12.5.csv` for J_obs = 12.5.
"""
waveform_filename(J_obs::Real) = @sprintf("phi_J%g.csv", J_obs)

"""
    save_waveforms(data, directory) -> Vector{String}

Write one file per fixed-areal-radius observer, named by `waveform_filename`,
with header `t,phi` and one row per stored sample: coordinate time t and
φ(t, J_obs). J_obs is the absolute areal radius from `Config.observers`
(in units of M0 when M0 = 1). The same samples are in observers.csv, which also
has the observer time T = t - f. Returns the paths written.
"""
function save_waveforms(data::RunData, directory::AbstractString)
    names = [waveform_filename(J) for J in data.settings.observers]
    # %g keeps 6 significant digits, so nearly equal radii could share a name.
    allunique(names) ||
        error("Observer radii $(data.settings.observers) give clashing waveform file names $names")
    mkpath(directory)
    paths = String[]
    for (j, name) in enumerate(names)
        path = joinpath(directory, name)
        open(path, "w") do io
            println(io, "t,phi")
            for i in eachindex(data.times)
                @printf(io, "%.17g,%.17g\n", data.times[i], data.signals[j][i])
            end
        end
        push!(paths, path)
    end
    return paths
end

"""Write portable, versioned CSV+TOML run data. Snapshot rows contain
(t,f,r,J,Jphi); fields are in geometrized units with M0 from config.toml.
Also writes the per-observer (t, phi) files of `save_waveforms`, replacing any
phi_J<J>.csv files of other observers in `directory`; `load_run` does not need
them. No GLMakie dependency and no global mutable output path.
"""
function save_run(data::RunData, directory::AbstractString)
    mkpath(directory)
    config=Dict(string(k)=>getfield(data.settings,k) for k in fieldnames(Config))
    open(joinpath(directory,"config.toml"),"w") do io
        TOML.print(io,config)
    end
    open(joinpath(directory,"metadata.toml"),"w") do io
        TOML.print(io,Dict(
            "schema_version"=>SCHEMA_VERSION,
            "code"=>"MMIEF v0.1.0",
            "units"=>"G=c=1; M0 in config.toml",
            "time_coordinate"=>"t; T_asymptotic=t-f",
            "grid_coordinate"=>"r; areal radius J=r+f",
            # Describes the baseline solver; scripts/run_constrained.jl
            # overwrites this label for constraint-solved runs.
            "evolution"=>"partially constrained, RK4+Hamiltonian projection",
            "snapshot_quantity"=>"Jphi=J*phi",
            "column_order"=>"a,B,Phi,Pi,phi"
        ))
    end
    open(joinpath(directory,"diagnostics.csv"),"w") do io
        println(io,"t,T_asymptotic,f,M_h,M_outer,hamiltonian_L2,scalar_definition_L2,horizon_error,horizon_preprojection_error")
        for i in eachindex(data.times)
            fields=(data.times[i],data.times[i]-data.shifts[i],data.shifts[i],
                    data.masses_h[i],data.masses_out[i],data.ham_l2[i],
                    data.scalar_constraint_l2[i],data.horizon_error[i],
                    data.horizon_preprojection[i])
            println(io,join((@sprintf("%.17g",x) for x in fields),","))
        end
    end
    # Long format: one row per (sample, observer) pair.
    open(joinpath(directory,"observers.csv"),"w") do io
        println(io,"t,T_asymptotic,J_obs,phi")
        for i in eachindex(data.times), j in eachindex(data.signals)
            fields=(data.times[i], data.times[i]-data.shifts[i],
                    data.settings.observers[j], data.signals[j][i])
            println(io,join((@sprintf("%.17g",x) for x in fields),","))
        end
    end
    # One row per (frame, retained grid point); J = r + f at that frame.
    open(joinpath(directory,"snapshots.csv"),"w") do io
        println(io,"t,f,r,J,Jphi")
        grid=view(data.r,1:data.settings.profile_stride:length(data.r))
        for k in eachindex(data.frame_times)
            profile=data.frame_profiles[k]
            length(profile)==length(grid) || error("Snapshot/grid length mismatch")
            t=data.frame_times[k]; f=data.frame_shift[k]
            for i in eachindex(profile)
                fields=(t,f,grid[i],grid[i]+f,profile[i])
                println(io,join((@sprintf("%.17g",x) for x in fields),","))
            end
        end
    end
    # Final fields on the full grid. f is not stored here; take it from the
    # last row of diagnostics.csv.
    open(joinpath(directory,"final_state.csv"),"w") do io
        println(io,"r,a,B,Phi,Pi,phi")
        for i in eachindex(data.r)
            fields=(data.r[i],data.final_state[i,1],data.final_state[i,2],
                    data.final_state[i,3],data.final_state[i,4],data.final_state[i,5])
            println(io,join((@sprintf("%.17g",x) for x in fields),","))
        end
    end
    # Convenience copies of the observer series: phi_J<J>.csv with (t, phi).
    # The directory describes one run, so first remove phi_J<J>.csv files left
    # by an earlier save_run with other observers (the files above are
    # overwritten anyway); no other file names are touched.
    current=Set(waveform_filename(J) for J in data.settings.observers)
    for name in readdir(directory)
        if occursin(r"^phi_J[-+.0-9eE]+\.csv$",name) && !(name in current)
            rm(joinpath(directory,name))
        end
    end
    save_waveforms(data,directory)
    return abspath(directory)
end

"""Read the versioned output independent of the simulation process.
Uses only Julia stdlib; no GLMakie or simulation required to load the data.
"""
function load_run(directory::AbstractString)
    meta=TOML.parsefile(joinpath(directory,"metadata.toml"))
    get(meta,"schema_version",nothing)==SCHEMA_VERSION ||
        error("Unsupported MMIEF output schema version")
    # Rebuild the Config from the saved keyword values.
    raw=TOML.parsefile(joinpath(directory,"config.toml"))
    opts=(Symbol(k)=>v for (k,v) in raw)
    p=Config(;opts...)
    # Grid and final state.
    data_rows=readlines(joinpath(directory,"final_state.csv"))
    length(data_rows)>=2 || error("No final state stored")
    r=Float64[]
    u=Matrix{Float64}(undef,length(data_rows)-1,5)
    for (i,line) in enumerate(data_rows[2:end])
        values=parse.(Float64,split(line,','))
        length(values)==6 || error("Invalid final_state.csv at row $i")
        push!(r,values[1])
        u[i,:].=values[2:6]
    end
    data=RunData(p,r,u)
    # Time series (header line skipped).
    for line in Iterators.drop(eachline(joinpath(directory,"diagnostics.csv")),1)
        v=parse.(Float64,split(line,','))
        length(v)==9 || error("Invalid diagnostics.csv")
        push!(data.times,v[1]); push!(data.shifts,v[3]); push!(data.masses_h,v[4])
        push!(data.masses_out,v[5]); push!(data.ham_l2,v[6])
        push!(data.scalar_constraint_l2,v[7]); push!(data.horizon_error,v[8])
        push!(data.horizon_preprojection,v[9])
    end
    # Observer rows are matched to Config.observers by their J_obs value.
    rows=0
    for line in Iterators.drop(eachline(joinpath(directory,"observers.csv")),1)
        v=parse.(Float64,split(line,','))
        length(v)==4 || error("Invalid observers.csv")
        j=findfirst(x->isapprox(x,v[3];atol=1e-10,rtol=0),p.observers)
        j===nothing && error("Unknown observer J=$(v[3])")
        push!(data.signals[j],v[4]); rows+=1
    end
    rows==length(data.times)*length(p.observers) || error("Missing observer rows")
    # Snapshots: a new frame starts whenever the time column changes.
    for line in Iterators.drop(eachline(joinpath(directory,"snapshots.csv")),1)
        v=parse.(Float64,split(line,','))
        length(v)==5 || error("Invalid snapshots.csv")
        if isempty(data.frame_times) || v[1]!=data.frame_times[end]
            push!(data.frame_times,v[1]); push!(data.frame_shift,v[2])
            push!(data.frame_profiles,Float32[])
        end
        push!(data.frame_profiles[end],Float32(v[5]))
    end
    expected=length(1:p.profile_stride:length(r))
    all(length(profile)==expected for profile in data.frame_profiles) ||
        error("Snapshot profile length mismatch")
    data.final_state=copy(u)
    return data
end
end # module DataOutput
