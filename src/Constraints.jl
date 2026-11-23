# Einstein constraints in MMIEF gauge (3+1, G = 1, ∂_r J = 1):
#   Hamiltonian: ∂_r a = -(a³-a)/(2J) - (a³J/2) B(2A+B) + 2πJa(Φ²+Π²)
#   Momentum:    ∂_r B = (A-B)/J + 4πΦΠ/a
# The momentum constraint is solved algebraically for A at every RHS call;
# the Hamiltonian constraint is enforced through the Misner–Sharp mass.
module Constraints
using ..Types: Config, IA, IB, IPHI, IPI
export mass_rhs, project_hamiltonian!, reconstruct_A!, hamiltonian_residual
"""
Momentum constraint reconstruction of radial extrinsic curvature.

A = B + J ∂_r B - 4πJΦΠ/a. Because A is built from the same discrete ∂_r B
(`dB`) used elsewhere, the discrete momentum residual is zero by construction.
"""
function reconstruct_A!(A,u,J,dB)
    @inbounds for i in eachindex(J)
        A[i] = u[i,IB] + J[i]*dB[i] -
               4π*J[i]*u[i,IPHI]*u[i,IPI]/u[i,IA]
    end
    return A
end
"""
Independent Hamiltonian residual from the Einstein constraint.

Diagnostic only: `ap` is a finite-difference ∂_r a, whereas the projection
below uses the radial mass ODE, so this residual is not zero by construction.
"""
@inline function hamiltonian_residual(J,a,ap,B,A,Phi,Pi)
    return ap+(a^3-a)/(2J)+(a^3*J/2)*B*(2*A+B)-2π*J*a*(Phi^2+Pi^2)
end
"""Radial Misner--Sharp mass equation for constraint projection.
This uses J'=1; 1/a^2=1-2m/J+(JB)^2.

∂_r m = 4πJ²[(Φ²+Π²)/(2a²) + (JB/a)ΦΠ], with a eliminated through
1/a² = 1 - 2m/J + (JB)² (positive root a > 0). Combining the definition of m
with the Hamiltonian and momentum constraints gives exactly this equation.
"""
@inline function mass_rhs(J::Float64, m::Float64, B::Float64,
                          Phi::Float64, Pi::Float64)
    ainvsq = 1 - 2m/J + (J*B)^2
    (isfinite(ainvsq) && ainvsq > 0) || error("Hamiltonian projection: nonpositive 1/a² at J=$J")
    ainv = sqrt(ainvsq)
    return 4π*J^2 * (0.5*(Phi^2+Pi^2)*ainvsq + J*B*ainv*Phi*Pi)
end

"""
Project the Hamiltonian constraint on each RK4 stage by integrating the
Misner--Sharp mass outward from m(J_h)=J_h/2, then reconstructing a.
4th-order radial RK with cubic midpoint interpolation, in vacuum exact.
This is a constrained-RK4 modification of the paper's original method.

The inner value m(J_h) = J_h/2 is the apparent-horizon condition aJB = 1
(for B_h > 0). Only column `a` is overwritten; B, Φ, Π are inputs.
"""
function project_hamiltonian!(u::Matrix{Float64}, f::Float64,
                              p::Config, r::Vector{Float64})
    n=length(r)
    m=(r[1]+f)/2       # horizon: m = J_h/2
    h=p.dr             # radial step; dJ = dr because J = r + f
    @inbounds for i in 1:n
        J=r[i]+f
        if i > 1
            # One radial RK4 step for m from J0 = J_{i-1} to J = J_i.
            j=i-1
            J0=r[j]+f
            B0, P0, Pi0 = u[j,IB], u[j,IPHI], u[j,IPI]
            B1, P1, Pi1 = u[i,IB], u[i,IPHI], u[i,IPI]
            if 2 <= j <= n-2
                # Fourth-order midpoint interpolation of stage fields.
                Bmid=(-u[j-1,IB]+9B0+9B1-u[j+2,IB])/16
                Pmid=(-u[j-1,IPHI]+9P0+9P1-u[j+2,IPHI])/16
                Pimid=(-u[j-1,IPI]+9Pi0+9Pi1-u[j+2,IPI])/16
            else
                # First and last interval: two-point average (second order).
                Bmid=(B0+B1)/2
                Pmid=(P0+P1)/2
                Pimid=(Pi0+Pi1)/2
            end
            k1=mass_rhs(J0,m,B0,P0,Pi0)
            k2=mass_rhs(J0+h/2,m+h*k1/2,Bmid,Pmid,Pimid)
            k3=mass_rhs(J0+h/2,m+h*k2/2,Bmid,Pmid,Pimid)
            k4=mass_rhs(J,m+h*k3,B1,P1,Pi1)
            # `2k2` means 2*k2 (k is not an exponent letter, unlike `2f2`).
            m+=h*(k1+2k2+2k3+k4)/6
        end
        # Reconstruct a > 0 from the mass and the evolved B.
        inva2=1-2m/J+(J*u[i,IB])^2
        (isfinite(inva2) && inva2 > 0) || error("Hamiltonian projection failed at J=$J")
        u[i,IA]=inv(sqrt(inva2))
    end
    return nothing
end

end # module Constraints
