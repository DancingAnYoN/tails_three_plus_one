# Massless Klein–Gordon field: initial pulse profile and the semi-discrete
# right-hand side of the first-order system for (Φ, Π, φ).
module ScalarField
using ..Types: Config, IPHI, IPI, IFIELD
export gaussian, scalar_rhs!
"""
Gaussian pulse multiplied by a C-infinity cutoff; return (phi, dphi/dr).

φ(r) = ε χ(r) exp[-(r - r0)²/σ²], where χ rises smoothly from 0 at
`cutoff_start` to 1 at `cutoff_end`, so φ vanishes identically near the horizon.
At t = 0, f = 0 and the grid coordinate r equals the areal radius J.
"""
function gaussian(r::Float64, p::Config)
    q = (r - p.pulse_r0) / p.pulse_sigma
    g = p.pulse_amplitude * exp(-q*q)
    gp = -2*q/p.pulse_sigma * g
    if r <= p.cutoff_start
        return 0.0, 0.0
    elseif r >= p.cutoff_end
        return g, gp
    end
    w = p.cutoff_end-p.cutoff_start
    x = (r-p.cutoff_start)/w
    # Standard smooth bump: chi=exp(-1/x)/(exp(-1/x)+exp(-1/(1-x))).
    e1 = exp(-1/x)
    e2 = exp(-1/(1-x))
    chi = e1 / (e1+e2)
    # dχ/dr = χ(1-χ)(1/x² + 1/(1-x)²)/w, by the chain rule through x(r).
    chip = chi*(1-chi)*(1/x^2 + 1/(1-x)^2)/w
    return chi*g, chi*gp+chip*g
end

"""
Scalar-field portion of the semi-discrete wave equations.
Flux arrays and their derivatives are precomputed by the coupled RHS.

With α/a = 1 - β substituted (MMIEF gauge):

    ∂_t Φ = ∂_r[βΦ + (1-β)Π]                          (dFFlux)
    ∂_t Π = J⁻² ∂_r[J²(βΠ + (1-β)Φ)] - 2ḟ Π/J          (dPFlux)
    ∂_t φ = βΦ + (1-β)Π                               (fFlux)

The Π equation is the conservative form ∂_t(J²Π) = ∂_r[J²(βΠ + (1-β)Φ)]
expanded with ∂_t J = ḟ (`fd`).
"""
function scalar_rhs!(du, u, J, fd, dFFlux, dPFlux, fFlux)
    @inbounds for i in eachindex(J)
        du[i,IPHI] = dFFlux[i]
        du[i,IPI] = dPFlux[i]/J[i]^2 - 2*fd*u[i,IPI]/J[i]
        du[i,IFIELD] = fFlux[i]
    end
    return nothing
end
end # module ScalarField
