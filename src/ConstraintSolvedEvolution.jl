module ConstraintSolvedEvolution

# Constraint-enforced MMIEF scalar-field evolution.
#
# The radial metric coefficient a is reconstructed from the
# Hamiltonian constraint at every RK4 stage.
# The momentum constraint determines A from B'.
# B remains an evolved Einstein variable.
#
# Evolved state: B, Φ, Π, φ (columns of u) and the scalar f. Column a is kept
# in u for compatibility with the shared modules, but its RHS is always zero
# and its value is overwritten by project_constraints!. This is constraint-
# solved, not "fully constrained" (see CONSTRAINED_FORMULATION.md §1).

using Printf
using ..Types: Config, IA, IB, IPHI, IPI, IFIELD
using ..Grid: radial_grid
using ..FiniteDifferences: deriv4!, ko6!
using ..Geometry: schwarzschild_fields, fdot_horizon
using ..Constraints: project_hamiltonian!, reconstruct_A!
using ..ScalarField: scalar_rhs!
using ..EinsteinEquations: gravity_rhs!
using ..InitialData: initial_data
using ..BoundaryConditions: outer_projection!, outer_rhs!
using ..Integrators: RK4Workspace, rk4_step!
using ..MethodOfLines: cfl_dt
using ..Observers: sample_fixed_J
using ..Diagnostics: diagnostics
using ..DataOutput: RunData, store_sample!, store_frame!

export run_constraint_solved, project_constraints!, constrained_rhs!,
       ConstrainedWorkspace, vacuum_constrained_rhs

# Preallocated RHS arrays; same meaning as Evolution.Workspace ("d" = ∂_r).
# Only the B equation is well-balanced here, because a is not evolved.
mutable struct ConstrainedWorkspace
    r::Vector{Float64}          # grid coordinates
    J::Vector{Float64}          # areal radius r + f
    A::Vector{Float64}          # K^r_r from the momentum constraint
    beta::Vector{Float64}       # shift β
    da::Vector{Float64}         # ∂_r a (only fed to outer_rhs!, whose a row is then zeroed)
    dB::Vector{Float64}         # ∂_r B
    dPhi::Vector{Float64}       # ∂_r Φ
    dBeta::Vector{Float64}      # ∂_r β
    dABeta::Vector{Float64}     # ∂_r(aβ)
    dFFlux::Vector{Float64}     # ∂_r[βΦ + (1-β)Π]
    dPFlux::Vector{Float64}     # ∂_r[J²(βΠ + (1-β)Φ)]
    aBeta::Vector{Float64}      # aβ
    fFlux::Vector{Float64}      # βΦ + (1-β)Π  (= ∂_t φ)
    pFlux::Vector{Float64}      # J²(βΠ + (1-β)Φ)
    vacuum_B_rhs::Vector{Float64}  # discrete Schwarzschild ∂_t B (well-balancing)
end

function ConstrainedWorkspace(r::Vector{Float64})
    arrays = [zeros(length(r)) for _ in 1:14]
    return ConstrainedWorkspace(r, arrays...)
end

"""Solve the Hamiltonian and momentum constraints at a time slice/RK4 stage.

The spatial Hamiltonian ODE is solved via the Misner--Sharp mass, starting
at the exact apparent-horizon boundary m(r_h)=J_h/2. Momentum constraint
reconstructs A pointwise when RHS is evaluated. The incoming outer scalar
characteristic is iterated with the constraint solve to remove the local
boundary inconsistency present in a single pass.

This operation changes `a` and the *outermost* `Pi`, but does not change
`B`, `Phi`, `phi`, or the horizon shift `f`.
"""
function project_constraints!(u::Matrix{Float64}, f::Float64,
                              p::Config, r::Vector{Float64})
    J_h = r[1] + f
    (isfinite(J_h) && J_h > 0 && u[1,IB] > 0) ||
        error("Invalid horizon: J_h must be positive and B_h must be positive")
    last_a = NaN
    last_pi = NaN
    converged = false
    # The two projections are coupled: the outer Π enters the mass integral
    # (and hence a at r_max), and a at r_max enters the outer scalar condition.
    # Alternate them until the outer values stop changing.
    for _ in 1:12
        project_hamiltonian!(u, f, p, r)
        outer_projection!(u, f, p, r)
        a = u[end,IA]
        pi = u[end,IPI]
        if isfinite(last_a) &&
           abs(a-last_a) <= 2e-12*max(1.0,abs(a)) &&
           abs(pi-last_pi) <= 2e-12*max(1.0,abs(pi))
            # The next Hamiltonian pass changes only floating point noise.
            project_hamiltonian!(u, f, p, r)
            converged = true
            break
        end
        last_a, last_pi = a, pi
    end
    converged || error("Outer Sommerfeld/Hamiltonian projection did not converge")
    horizon_error = u[1,IA]*J_h*u[1,IB] - 1
    abs(horizon_error) < 1e-10 ||
        error("Horizon mass boundary and a*J*B=1 inconsistent: $horizon_error")
    return nothing
end

"""Semi-discrete RHS for (B, Phi, Pi, phi, f).

Column `a` has identically zero RHS; its value is freshly solved from the
Hamiltonian constraint before each evaluation. `A` is obtained from the
momentum constraint, not evolved. `f` uses horizon locking.

The subtracted stationary Schwarzschild B residual is a discrete
well-balancing correction, *not* an additional field equation.
"""
function constrained_rhs!(du::Matrix{Float64},u::Matrix{Float64},f::Float64,
                          p::Config,w::ConstrainedWorkspace;
                          well_balance::Bool=true)
    n = length(w.r)
    @views begin
        a = u[:,IA]
        B = u[:,IB]
        Phi = u[:,IPHI]
        Pi = u[:,IPI]
        deriv4!(w.da,a,p.dr)
        deriv4!(w.dB,B,p.dr)
        deriv4!(w.dPhi,Phi,p.dr)
        fd=fdot_horizon(u,f,p)
        (isfinite(fd) && 0 <= fd < 1) ||
            error("MMIEF horizon locking failed: f_dot=$fd (requires 0 <= f_dot < 1)")
        @inbounds for i in 1:n
            Ji = w.r[i]+f
            ai = a[i]
            (isfinite(Ji) && Ji>0 && isfinite(ai) && ai>0) ||
                error("Invalid geometry at grid index $i")
            w.J[i] = Ji
        end
        # A from the momentum constraint; then shift and scalar fluxes
        # (the lapse is never formed: α/a = 1-β).
        reconstruct_A!(w.A,u,w.J,w.dB)
        @inbounds for i in 1:n
            Ji=w.J[i]
            q=a[i]*Ji*B[i]
            den=1+q
            den>0 || error("MMIEF gauge denominator 1+aJB <=0 at index $i")
            b=(fd+q)/den
            w.beta[i]=b
            w.aBeta[i]=a[i]*b
            w.fFlux[i]=b*Phi[i]+(1-b)*Pi[i]
            w.pFlux[i]=Ji^2*(b*Pi[i]+(1-b)*Phi[i])
        end
        deriv4!(w.dBeta,w.beta,p.dr)
        deriv4!(w.dABeta,w.aBeta,p.dr)
        deriv4!(w.dFFlux,w.fFlux,p.dr)
        deriv4!(w.dPFlux,w.pFlux,p.dr)
        # Interior equations, KO dissipation, then the outer-boundary row.
        gravity_rhs!(du,u,w.J,w.A,w.beta,w.dB,w.dBeta,w.dABeta)
        scalar_rhs!(du,u,w.J,fd,w.dFFlux,w.dPFlux,w.fFlux)
        ko6!(du,u,p.dr,p.ko_strength)
        outer_rhs!(du,u,w.J,w.beta,fd,w.da,w.dB,w.dPhi,w.dBeta)
        # `a` is *always* constrained: no independent gravitational update.
        du[:,IA] .= 0.0
        if well_balance
            @inbounds for i in 1:n
                du[i,IB] -= w.vacuum_B_rhs[i]
            end
        end
        return fd
    end
end

"""Construct and cancel the stationary Schwarzschild discrete B residual."""
function set_vacuum_balance!(p::Config,w::ConstrainedWorkspace)
    N=length(w.r)
    u=zeros(N,5)
    du=similar(u)
    for i in eachindex(w.r)
        s=schwarzschild_fields(w.r[i],p.M0)
        u[i,IA]=s.a
        u[i,IB]=s.B
    end
    project_constraints!(u,0.0,p,w.r)
    constrained_rhs!(du,u,0.0,p,w;well_balance=false)
    @views w.vacuum_B_rhs .= du[:,IB]
    return nothing
end

"""
Vacuum check for the constraint-solved RHS: returns `(maximum_rhs, fdot,
horizon_error)` for exact Schwarzschild data. Well-balancing is on, so a
near-zero `maximum_rhs` checks the implementation, not the truncation error.
"""
function vacuum_constrained_rhs(;dr::Float64=0.5, rmax::Float64=70.0)
    p=Config(dr=dr,rmax=rmax,tmax=1.0,pulse_amplitude=0.0)
    r=radial_grid(p)
    u=zeros(length(r),5)
    for i in eachindex(r)
        s=schwarzschild_fields(r[i],p.M0)
        u[i,IA]=s.a
        u[i,IB]=s.B
    end
    w=ConstrainedWorkspace(r)
    set_vacuum_balance!(p,w)
    project_constraints!(u,0.0,p,r)
    du=similar(u)
    fd=constrained_rhs!(du,u,0.0,p,w)
    return (maximum_rhs=maximum(abs,du),fdot=fd,
            horizon_error=u[1,IA]*r[1]*u[1,IB]-1)
end

"""Evolve a constrained 3+1 Einstein--Klein--Gordon MMIEF model.

The freely advanced radial fields are B, Phi, Pi, phi; the scalar ODE is f.
The `a` column is reconstructed on each of four RK4 stages. `A` is obtained
from the momentum constraint at RHS evaluation.

Returns the existing `RunData` type, compatible with all `viz/` scripts.
"""
function run_constraint_solved(p::Config;verbose::Bool=true)
    r,u=initial_data(p)
    f=0.0
    project_constraints!(u,f,p,r)
    w=ConstrainedWorkspace(r)
    set_vacuum_balance!(p,w)
    rk=RK4Workspace(u)
    eval_rhs! = (du,x,shift)->constrained_rhs!(du,x,shift,p,w)
    # project_constraints! already includes the outer scalar projection, so
    # the separate boundary hook is a no-op.
    project! = (x,shift)->project_constraints!(x,shift,p,r)
    noop! = (x,shift)->nothing
    prior_residual = (x,shift)->x[1,IA]*(r[1]+shift)*x[1,IB]-1
    diag0=diagnostics(u,f,r,p)
    outer_mass0=diag0.M_out
    # Fail-fast tolerance for unphysical horizon growth (not an error estimate).
    budget_tol=max(p.mass_budget_fraction*p.M0,
                   0.1*max(0.0,outer_mass0-p.M0))
    data=RunData(p,r,u)
    store_sample!(data,0.0,u,f)
    store_frame!(data,0.0,u,f)
    if verbose
        @printf("Constrained MMIEF: N=%d dr=%g tmax=%g | M_out=%1.9g H_L2=%.3e\n",
                length(r),p.dr,p.tmax,outer_mass0,diag0.ham_l2)
    end
    t=0.0; step=0
    next_sample=p.output_dt
    next_frame=p.snapshot_dt
    next_report=25*p.M0
    while t < p.tmax-1e-12
        # Time step from the largest coordinate speed: scalar speeds -1 and
        # 1-2β, plus the shift β that advects B.
        fd=fdot_horizon(u,f,p)
        vmax=1.0
        @inbounds for i in eachindex(r)
            q=u[i,IA]*(r[i]+f)*u[i,IB]
            den=1+q
            den>0 || error("Gauge denominator failed while computing CFL")
            beta=(fd+q)/den
            vmax=max(vmax,abs(1-2*beta),abs(beta))
        end
        dt=cfl_dt(p.dr,p.cfl,vmax,p.tmax-t)
        f,preproj=rk4_step!(u,f,dt,rk;rhs! = eval_rhs!,
                            project! = project!, boundary! = noop!,
                            preprojection=prior_residual)
        t+=dt; step+=1
        if !isfinite(f) || any(x->!isfinite(x),u)
            error("Nonfinite constrained solution at time $t")
        end
        if p.check_mass_budget && (r[1]+f)/2 > outer_mass0+budget_tol
            error(@sprintf("Horizon mass budget exceeded at t=%.6g: Mh=%.6g, initial outer mass=%.6g",
                           t,(r[1]+f)/2,outer_mass0))
        end
        for Jobs in p.observers
            # This interpolator fails helpfully if the observer crosses the horizon.
            if !(r[1] <= Jobs-f <= r[end])
                sample_fixed_J(u,r,f,Jobs,p;time=t)
            end
        end
        if t+1e-10>=next_sample || t>=p.tmax-1e-12
            store_sample!(data,t,u,f,preproj)
            next_sample+=p.output_dt
        end
        if t+1e-10>=next_frame || t>=p.tmax-1e-12
            store_frame!(data,t,u,f)
            next_frame+=p.snapshot_dt
        end
        if verbose && (t+1e-10>=next_report || t>=p.tmax-1e-12)
            d=diagnostics(u,f,r,p)
            @printf("t/M0=%7.2f Mh=%1.9g f=%9.3e H_L2=%9.2e horizon=%9.2e\n",
                    t/p.M0,d.M_h,f,d.ham_l2,d.horizon_error)
            next_report+=25*p.M0
        end
    end
    data.final_state=copy(u)
    return data
end

end # module ConstraintSolvedEvolution
