using Test
using MMIEF

@testset "Constraint-solved MMIEF: vacuum and short scalar test" begin
    v=MMIEF.vacuum_constrained_rhs(dr=0.5,rmax=70.0)
    @test v.maximum_rhs < 1e-11
    @test v.fdot == 0.0
    @test abs(v.horizon_error) < 1e-11

    vac=MMIEF.Config(rmax=70.0,dr=0.5,cfl=0.1,tmax=0.25,
        pulse_amplitude=0.0,output_dt=0.125,snapshot_dt=0.125,
        ko_strength=0.025)
    data=MMIEF.run_constraint_solved(vac;verbose=false)
    @test length(data.times)>=2
    @test maximum(abs.(data.shifts)) < 1e-11
    @test maximum(abs.(data.masses_h .- vac.M0)) < 1e-11
    @test maximum(abs.(data.horizon_error)) < 1e-10
    @test maximum(abs.(data.masses_out .- vac.M0)) < 1e-7
    @test all(isfinite, data.final_state)

    pulse=MMIEF.Config(rmax=70.0,dr=0.5,cfl=0.1,tmax=0.25,
        pulse_amplitude=1e-5,output_dt=0.125,snapshot_dt=0.125,
        ko_strength=0.025)
    d=MMIEF.run_constraint_solved(pulse;verbose=false)
    @test all(isfinite,d.final_state)
    @test all(isfinite,d.ham_l2)
    @test maximum(abs.(d.horizon_error)) < 1e-10
    @test maximum(d.ham_l2) < 1e-2
    @test d.masses_h[end] >= pulse.M0-1e-10
    @test !isempty(d.frame_times)
end
