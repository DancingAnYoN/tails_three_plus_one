# Approximate outgoing-wave conditions at the outer boundary r = r_max.
# There is no explicit inner boundary condition: the horizon r_h is an outflow
# boundary (scalar speeds -1 and 1-2β_h = -ḟ ≤ 0), handled by the one-sided
# stencils in FiniteDifferences.deriv4!. All conditions here are approximate at
# finite radius; reflections must be checked by increasing r_max.
module BoundaryConditions
using ..Types: Config, IA, IB, IPHI, IPI, IFIELD
using ..Geometry: fdot_horizon
export outer_projection!, outer_rhs!
"""
Project only the incoming scalar characteristic at the distant outer boundary.
Uses the asymptotic outgoing spherical-wave condition, Eq. (4.11).

At r_max the incoming characteristic variable is Π + Φ (speed -1). Its value
follows from the Sommerfeld condition (∂_t + c ∂_r)(Jφ) = 0 with outgoing speed
c = 1 - 2β, which gives (1-β)(Π + Φ) + (1 - 2β + ḟ) φ/J = 0. Only Π at the last
grid point is changed.
"""
function outer_projection!(u::Matrix{Float64}, f::Float64, p::Config, r::Vector{Float64})
    i = length(r)
    fd = fdot_horizon(u,f,p)
    q = u[i,IA]*(r[i]+f)*u[i,IB]
    den = 1+q
    den > 0 || error("Outer gauge denominator <= 0")
    beta = (fd+q)/den
    omb = 1-beta
    omb > 0 || error("Outer lapse factor <= 0")
    c = 1-2beta     # outgoing radial null speed
    u[i,IPI] = -u[i,IPHI] - (fd+c)*u[i,IFIELD]/((r[i]+f)*omb)
    return nothing
end

"""Outer radial RHS following the original code's approximate outgoing conditions.
Precomputed derivatives are passed explicitly; this routine is not the source
of their discretization.

Overwrites the RHS at the last grid point with outgoing advection conditions
(c = 1 - 2β), each a Sommerfeld condition on a quantity with the expected falloff:
  a: (∂_t + c∂_r)[J(a-1)] = 0     since a - 1 = O(1/J)
  B: (∂_t + c∂_r)[J²B]    = 0     since B = O(1/J²)
  Φ: radial derivative of the condition on Jφ, using ∂_r c = -2∂_r β.
The φ row keeps its interior RHS.
"""
function outer_rhs!(du,u,J,beta,fd,dAmetric,dB,dPhi,dBeta)
    i=size(u,1)
    Ji=J[i]; b=beta[i]; c=1-2*b
    du[i,IA] = -c*dAmetric[i]-(fd+c)*(u[i,IA]-1)/Ji
    du[i,IB] = -c*dB[i]-2*(fd+c)*u[i,IB]/Ji
    du[i,IPHI] = -c*dPhi[i]-(fd+c-2*Ji*dBeta[i])*u[i,IPHI]/Ji +
                  (fd+c+2*Ji*dBeta[i])*u[i,IFIELD]/Ji^2
    # Placeholder that keeps Π + Φ fixed within a stage; outer_projection!
    # then imposes its actual value after each stage.
    du[i,IPI] = -du[i,IPHI] # projection enforces actual Pi value
    return nothing
end
end # module BoundaryConditions
