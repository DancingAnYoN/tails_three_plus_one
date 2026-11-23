# Method-of-lines layer: bundles a semi-discrete RHS with its projection and
# boundary hooks, and chooses the time step. Knows nothing about the physics.
module MethodOfLines
using ..Integrators: RK4Workspace, rk4_step!
export MOLSystem, mol_step!, cfl_dt

"""The semi-discrete spatial RHS and numerical constraints/boundary hooks.
No dependency on any particular PDE, gravitational geometry, or plotting.
"""
struct MOLSystem{R,P,B,D}
    rhs!::R            # (du, u, f) -> ḟ
    project!::P        # (u, f) -> nothing; constraint projection of a stage
    boundary!::B       # (u, f) -> nothing; boundary projection of a stage
    preprojection::D   # (u, f) -> Float64; residual recorded before projection
end

"""Evolve one time step of an arbitrary method-of-lines semidiscretization."""
function mol_step!(u::Matrix{Float64},f::Float64,dt::Float64,
                   work::RK4Workspace,model::MOLSystem)
    return rk4_step!(u,f,dt,work;
                     rhs! = model.rhs!, project! = model.project!,
                     boundary! = model.boundary!, preprojection=model.preprojection)
end

"""
Conservative CFL selection with final-step truncation.

Δt = C Δr / maxspeed, shortened if needed so the run ends exactly at t_max.
"""
function cfl_dt(dr::Real,C::Real,maxspeed::Real,remaining::Real)
    dr>0 && C>0 && maxspeed>0 && remaining>0 ||
        throw(ArgumentError("Invalid CFL inputs"))
    return min(C*dr/maxspeed,remaining)
end
end # module MethodOfLines
