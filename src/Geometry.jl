# Gauge reconstruction, exact Schwarzschild data, Misner–Sharp mass and the
# horizon-locking law. Formulas are specialised to 3+1 with J = r + f(t).
module Geometry
using ..Types: Config, IA, IB, IPHI, IPI
export schwarzschild_fields, schwarzschild_Bprime, fdot_horizon,
       lapse_shift, misner_sharp_mass, horizon_radius, mmief_outgoing_delay
"""
Vacuum Schwarzschild values for f=0, useful as an exact stationary benchmark.

In vacuum the MMIEF gauge reduces to ingoing Eddington–Finkelstein (Kerr–Schild)
slicing: a = √(1+2M/J), β = 2M/(J+2M), α = a(1-β) = √(J/(J+2M)).
At the horizon J = 2M one has B = 1/(2√2 M), so aJB = 1.
"""
function schwarzschild_fields(J::Float64, M::Float64)
    a = sqrt(1 + 2M/J)
    den = (J*(J+2M))^1.5
    B = 2M*(J+2M)/den
    A = -2M*(J+M)/den
    beta = 2M/(J+2M)
    return (a=a, B=B, A=A, beta=beta)
end

"""
Analytic radial derivative of the Schwarzschild B component.

Used by `InitialData` to integrate only the deviation δB = B - B_Schw.
"""
@inline schwarzschild_Bprime(r::Float64,M::Float64) =
    -2M*(2r+3M)/(r^(5/2)*(r+2M)^(3/2))

"""
Gauge reconstruction using J=r+f and incoming null slicing.

From ∂_t J = β ∂_r J - αJB with J = r + f (so ∂_t f = β - αJB) and the
ingoing-null condition α = a(1-β):

    β = (ḟ + aJB)/(1 + aJB),   α = a(1-β) = a(1-ḟ)/(1 + aJB).

The RHS loops inline the same β formula; this function is used by tests.
"""
function lapse_shift(a::Real,J::Real,B::Real,fdot::Real)
    q = a*J*B
    denominator = 1+q
    denominator>0 || error("Invalid gauge denominator (1+aJB)<=0")
    beta=(fdot+q)/denominator
    alpha=a*(1-beta)
    alpha>0 || error("Nonpositive lapse")
    return (alpha=alpha,beta=beta)
end
# Misner–Sharp mass m = (J/2)(1 - g^{μν}∂_μJ ∂_νJ) = (J/2)(1 - 1/a² + J²B²),
# using ∂_r J = 1. An apparent horizon (aJB = 1) has m = J/2.
@inline misner_sharp_mass(J,a,B) = J/2*(1-1/a^2+(J*B)^2)
# Areal radius of the horizon, J_h = r_h + f with r_h = 2*M0; M_h = J_h/2.
@inline horizon_radius(p::Config,f::Real)=2*p.M0+f
"""
Horizon-locked coordinate growth rate, Eq. (2.20).

ḟ = 4π J_h² (Φ + Π)²/a_h², evaluated at grid row 1 (r = r_h). This is the unique
ḟ for which ∂_t(aJB) = 0 at r_h, given the evolution equations and both
constraints, so the apparent horizon aJB = 1 stays on the grid point r_h while
its areal radius J_h = r_h + f grows. Φ + Π is the ingoing scalar characteristic
variable: only ingoing flux feeds the hole, and ḟ ≥ 0.
"""
@inline function fdot_horizon(u::Matrix{Float64}, f::Float64, p::Config)
    J = 2*p.M0+f
    return 4π*(J/u[1,IA])^2*(u[1,IPHI]+u[1,IPI])^2
end

"""
    mmief_outgoing_delay(J, J_ref, M) -> ΔT

Observer-time delay of a stationary outgoing null ray from areal radius J_ref
to J on a Schwarzschild background of mass M. There β = 2M/(J+2M), so outgoing
rays move with dJ/dT = 1 - 2β = (J-2M)/(J+2M), and

    ΔT = (J - J_ref) + 4M ln[(J - 2M)/(J_ref - 2M)].

This is an ingoing Eddington–Finkelstein time delay: a difference of
Schwarzschild tortoise radii (coefficient 2M) measures Schwarzschild time and
is smaller by 2M ln[...]. Negative for J < J_ref.
"""
function mmief_outgoing_delay(J::Real, J_ref::Real, M::Real)
    (J > 2M && J_ref > 2M) ||
        throw(ArgumentError("Both radii must lie outside the horizon radius 2M"))
    return (J - J_ref) + 4M*log((J - 2M)/(J_ref - 2M))
end

end # module Geometry
