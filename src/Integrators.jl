# Generic classical RK4 for a matrix state plus one scalar ODE (the horizon
# shift f). Physics enters only through the callbacks.
module Integrators
export RK4Workspace, rk4_step!

"""A generic RK4 integrator for an array state coupled to one scalar ODE.
The physics RHS is supplied as a callback `rhs!(du,u,f)::Float64`.
Optional projection and boundary callbacks modify stage states in place.
This integrator has no dependence on Einstein equations or GLMakie.
"""
struct RK4Workspace
    k1::Matrix{Float64}
    k2::Matrix{Float64}
    k3::Matrix{Float64}
    k4::Matrix{Float64}
    stage::Matrix{Float64}    # scratch state U_n + c dt k for stages 2-4
end
RK4Workspace(u::Matrix{Float64}) =
    RK4Workspace(similar(u),similar(u),similar(u),similar(u),similar(u))

"""
Advance `u` and `f` by one RK4 step of size `dt`; return `(f_new, residual)`.

`rhs!(du, u, f)` fills du and returns ḟ. Each intermediate stage state is
passed through `project!` and then `boundary!` before the RHS is evaluated, and
the accepted state is projected the same way:

    U⁽ˢ⁾ = P(Uₙ + c_s dt k_{s-1}),  k_s = F(U⁽ˢ⁾),
    Uₙ₊₁ = P(Uₙ + dt/6 (k₁ + 2k₂ + 2k₃ + k₄)).

Stage 1 uses `u` as given (already projected at the end of the previous step).
`residual = preprojection(u, f_new)` is evaluated on the raw update before the
final projection, so it measures how far one unprojected step drifts.
"""
function rk4_step!(u::Matrix{Float64}, f::Float64, dt::Float64,
                   w::RK4Workspace;
                   rhs!, project! = (state,shift)->nothing,
                   boundary! = (state,shift)->nothing,
                   preprojection = (state,shift)->0.0)
    k1,k2,k3,k4,tmp=w.k1,w.k2,w.k3,w.k4,w.stage
    # Stage 1.
    f1 = rhs!(k1,u,f)
    # Stage 2: half step with k1.
    @inbounds for i in eachindex(u)
        tmp[i] = u[i]+dt*k1[i]/2
    end
    f2arg=f+dt*f1/2
    project!(tmp,f2arg)
    boundary!(tmp,f2arg)
    f2=rhs!(k2,tmp,f2arg)
    # Stage 3: half step with k2.
    @inbounds for i in eachindex(u)
        tmp[i]=u[i]+dt*k2[i]/2
    end
    f3arg=f+dt*f2/2
    project!(tmp,f3arg)
    boundary!(tmp,f3arg)
    f3=rhs!(k3,tmp,f3arg)
    # Stage 4: full step with k3.
    @inbounds for i in eachindex(u)
        tmp[i]=u[i]+dt*k3[i]
    end
    f4arg=f+dt*f3
    project!(tmp,f4arg)
    boundary!(tmp,f4arg)
    f4=rhs!(k4,tmp,f4arg)
    # Accepted update with the classical weights 1, 2, 2, 1.
    @inbounds for i in eachindex(u)
        u[i] += dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6
    end
    # Explicit * is essential in Julia: `2f2` is a Float32 literal (200).
    fnew=f+dt*(f1+2*f2+2*f3+f4)/6
    residual=preprojection(u,fnew)
    project!(u,fnew)
    boundary!(u,fnew)
    return fnew,residual
end
end # module Integrators
