using Test
using MMIEF

@testset "Versioned output and round-trip" begin
    p=Config(rmax=80.0,dr=0.5,tmax=1.0,pulse_amplitude=0.0,
             output_dt=0.5,snapshot_dt=0.5,profile_stride=2)
    data=run_evolution(p;verbose=false)
    mktempdir() do dir
        output=save_run(data,joinpath(dir,"case"))
        @test isfile(joinpath(output,"metadata.toml"))
        @test isfile(joinpath(output,"config.toml"))
        @test isfile(joinpath(output,"diagnostics.csv"))
        @test isfile(joinpath(output,"observers.csv"))
        @test isfile(joinpath(output,"snapshots.csv"))
        recovered=load_run(output)
        @test recovered.times==data.times
        @test recovered.shifts==data.shifts
        @test recovered.masses_h==data.masses_h
        @test recovered.signals==data.signals
        @test recovered.frame_times==data.frame_times
        @test recovered.frame_profiles==data.frame_profiles
        @test recovered.final_state==data.final_state
        @test recovered.settings.observers==data.settings.observers
        @test isfile(save_csv(data,joinpath(dir,"legacy.csv")))
        @test all(isfile(joinpath(output,name))
                  for name in ("phi_J10.csv","phi_J30.csv","phi_J50.csv"))
    end
end

@testset "Per-observer (t, phi) waveform files" begin
    # Non-vacuum run with non-default observers, so phi is nonzero and the
    # file names exercise the %g formatting (12.5 -> "12.5", 40.0 -> "40").
    p=Config(rmax=80.0,dr=0.5,tmax=1.0,output_dt=0.25,snapshot_dt=0.5,
             profile_stride=2,observers=[12.5,40.0])
    data=run_evolution(p;verbose=false)
    @test any(!iszero,data.signals[1])
    mktempdir() do dir
        output=save_run(data,joinpath(dir,"case"))
        @test waveform_filename.(p.observers)==["phi_J12.5.csv","phi_J40.csv"]
        for (j,J) in enumerate(p.observers)
            lines=readlines(joinpath(output,waveform_filename(J)))
            @test lines[1]=="t,phi"
            rows=[parse.(Float64,split(line,',')) for line in lines[2:end]]
            @test all(length(v)==2 for v in rows)
            # %.17g round-trips Float64 exactly.
            @test [v[1] for v in rows]==data.times
            @test [v[2] for v in rows]==data.signals[j]
        end
        @test isapprox(data.times[end],p.tmax;atol=1e-12)
        @test basename.(save_waveforms(data,joinpath(dir,"again")))==
              ["phi_J12.5.csv","phi_J40.csv"]
        # Re-saving into the same directory drops waveform files of observers
        # that are no longer in the run, and leaves unrelated files alone.
        write(joinpath(output,"phi_J99.csv"),"t,phi\n")
        write(joinpath(output,"notes.csv"),"keep\n")
        save_run(data,output)
        @test !isfile(joinpath(output,"phi_J99.csv"))
        @test isfile(joinpath(output,"notes.csv"))
        @test isfile(joinpath(output,"phi_J12.5.csv"))
    end
    # Radii that coincide to 6 significant digits would overwrite each other.
    q=Config(rmax=80.0,dr=0.5,tmax=1.0,observers=[30.0,30.0000001])
    clash=RunData(q,radial_grid(q),zeros(length(radial_grid(q)),5))
    mktempdir() do dir
        @test_throws ErrorException save_waveforms(clash,dir)
    end
end
