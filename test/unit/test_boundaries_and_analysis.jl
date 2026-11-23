using Test
using MMIEF

@testset "Outer scalar boundary projection" begin
    p=Config(rmax=80.0,dr=0.5,tmax=1.0,pulse_amplitude=0.0)
    r,u=initial_data(p)
    u[end,3]=0.003
    u[end,5]=-0.002
    MMIEF.BoundaryConditions.outer_projection!(u,0.0,p,r)
    fd=MMIEF.Geometry.fdot_horizon(u,0.0,p)
    a=u[end,1]; B=u[end,2]; J=r[end]
    beta=(fd+a*J*B)/(1+a*J*B)
    c=1-2*beta
    residual=(1-beta)*(u[end,4]+u[end,3])+(fd+c)*u[end,5]/J
    @test abs(residual)<1e-14
end

@testset "Synthetic QNM and tail recovery" begin
    t=collect(0.0:0.5:125.0)
    gamma=0.094
    omega=0.111
    y=exp.(-gamma.*t).*cos.(omega.*t.+0.4)
    fit=fit_qnm(t,y;tlo=45.0,thi=100.0,mass=1.0)
    @test fit !== nothing
    if fit !== nothing
        @test abs(fit.omegaM-omega)<0.015
        @test abs(fit.gammaM-gamma)<0.015
    end
    t2=collect(100.0:1.0:500.0)
    y2=t2.^(-3)
    pidx=local_power_index(t2,y2;peak_halfwidth=0.0,
                            fit_halfwidth=30.0,tstart=120.0)
    kept=[x for x in pidx if isfinite(x)]
    @test length(kept)>50
    @test abs(kept[end]-3.0)<0.05
end

@testset "Local log-log slope of a tail" begin
    # Exact for a pure power law, including truncated windows at both ends.
    t=collect(0.5:0.5:600.0)
    y=-2.5e-6 .* t.^(-3)
    p=local_loglog_slope(t,y;halfwidth=25.0)
    @test all(isfinite,p)
    @test maximum(abs.(p .- 3.0))<1e-10
    # tstart: no estimate before it, and windows do not reach below it.
    p=local_loglog_slope(t,y;halfwidth=25.0,tstart=100.0)
    @test all(isnan,p[t .< 100.0]) && all(isfinite,p[t .>= 100.0])
    # Undefined in any window that contains a sign change of y.
    z=t.^(-3) .* (t .- 200.0)
    p=local_loglog_slope(t,z;halfwidth=10.0)
    @test all(isnan,p[abs.(t .- 200.0) .<= 10.0])
    @test all(isfinite,p[t .> 211.0])
    # The peak-envelope index is biased high by about t/(t-h) on this tail.
    q=local_power_index(t,abs.(y);peak_halfwidth=12.0,fit_halfwidth=24.0,tstart=100.0)
    i=findfirst(==(300.0),t)
    @test isapprox(q[i],3*300/(300-12);rtol=2e-3)
end

@testset "Tail model fit A t^-γ exp(b1/t + b2/t²)" begin
    t=collect(0.0:0.5:1000.0)
    model(A,g,b1,b2)=A .* t.^(-g) .* exp.(b1 ./ t .+ b2 ./ t.^2)
    # Exact recovery of all parameters, including a negative amplitude.
    y=model(-3.2e-4,3.0,95.0,1800.0)
    f=fit_tail_model(t,y;tmin=150.0)
    @test isapprox(f.gamma,3.0;atol=1e-8)
    @test isapprox(f.b1,95.0;rtol=1e-6) && isapprox(f.b2,1800.0;rtol=1e-5)
    @test isapprox(f.A,-3.2e-4;rtol=1e-7) && f.rms<1e-10
    @test f.n==count(x->150<=x<=1000,t) && f.tmin==150.0 && f.tmax==1000.0
    # Lower orders: a pure power law is recovered exactly with order 0 and 1.
    p=model(2.0,2.5,0.0,0.0)
    @test isapprox(fit_tail_model(t,p;tmin=150.0,order=0).gamma,2.5;atol=1e-10)
    f1=fit_tail_model(t,p;tmin=150.0,tmax=600.0,order=1)
    @test isapprox(f1.gamma,2.5;atol=1e-8) && abs(f1.b1)<1e-5 && f1.tmax==600.0
    # A sign change in the window is rejected.
    @test_throws ErrorException fit_tail_model(t,(t .- 300.0) .* t.^(-4);tmin=150.0)
end

@testset "MMIEF outgoing delay" begin
    # Direct Simpson integration of dT/dJ = (J+2M)/(J-2M) along an outgoing ray.
    M=1.3; J1=10.0; J2=50.0; n=20000; h=(J2-J1)/n
    g(J)=(J+2M)/(J-2M)
    integral=h/3*(g(J1)+g(J2)+4*sum(g(J1+(2k-1)*h) for k in 1:n÷2)+
                  2*sum(g(J1+2k*h) for k in 1:n÷2-1))
    @test isapprox(mmief_outgoing_delay(J2,J1,M),integral;rtol=1e-10)
    @test mmief_outgoing_delay(J1,J2,M) ≈ -integral
    # Larger than the Schwarzschild tortoise difference by 2M ln[(J2-2M)/(J1-2M)].
    rstar(J)=J+2M*log(J/(2M)-1)
    @test mmief_outgoing_delay(J2,J1,M)-(rstar(J2)-rstar(J1)) ≈ 2M*log((J2-2M)/(J1-2M))
    @test_throws ArgumentError mmief_outgoing_delay(2.0,30.0,1.0)
end

@testset "Joint ringdown + tail fit" begin
    # Fundamental scalar l=0 mode plus a tail template of comparable size late
    # in the window, as in a J = 30 run.
    t=collect(85.0:0.5:200.0); tau=t .- t[1]
    w,G=0.110455,0.104896
    mode=exp.(-G .* tau) .* (1.5e-6 .* cos.(w .* tau) .- 0.4e-6 .* sin.(w .* tau))
    template=0.054 .* t.^(-3) .* exp.(70.0 ./ t .+ 5900.0 ./ t.^2)
    y=mode .+ 0.8 .* template
    f=fit_qnm_tail(t,y;tlo=85.0,thi=200.0,tail=template)
    @test abs(f.omegaM-w)<1e-4 && abs(f.gammaM-G)<1e-4
    @test isapprox(f.kappa,0.8;rtol=1e-3) && isapprox(f.a,1.5e-6;rtol=1e-3)
    @test f.rel_rms<1e-3 && f.n==length(t)
    # Without the tail term the same data give a badly biased mode.
    f0=fit_qnm_tail(t,y;tlo=85.0,thi=200.0)
    @test isnan(f0.kappa) && abs(f0.gammaM/G-1)>0.05
    # Physical frequencies scale with the mass.
    fm=fit_qnm_tail(t .* 2,y;tlo=170.0,thi=400.0,tail=template,mass=2.0)
    @test abs(fm.omegaM-w)<1e-4 && abs(fm.gammaM-G)<1e-4
    @test fit_qnm_tail(t,y;tlo=85.0,thi=88.0) === nothing
end
