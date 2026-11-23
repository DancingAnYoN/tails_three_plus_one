# Evolution equations for the gravitational fields a and B in MMIEF gauge.
# A = K^r_r is not evolved: the caller reconstructs it from the momentum
# constraint (Constraints.reconstruct_A!) and passes it in.
module EinsteinEquations
using ..Types: IA, IB
export gravity_rhs!
"""
The two evolved Einstein equations in the MMIEF gauge.

    ∂_t a = -a²(1-β)A + ∂_r(aβ)
    ∂_t B = β ∂_r B + a(1-β)B(A + 2B) + (1-β)(a - 1/a)/J² + ∂_r β/(aJ)

These are the ADM equations with J = r + f and α = a(1-β) substituted; all
derivatives of the lapse have been eliminated. Inputs: `dB` = ∂_r B,
`dBeta` = ∂_r β, `dABeta` = ∂_r(aβ). The constraint-solved path zeroes the
a row afterwards because there a is reconstructed, not evolved.
"""
function gravity_rhs!(du,u,J,A,beta,dB,dBeta,dABeta)
    @inbounds for i in eachindex(J)
        ai=u[i,IA]; Bi=u[i,IB]; Ai=A[i]; b=beta[i]; Ji=J[i]
        du[i,IA] = -ai^2*(1-b)*Ai + dABeta[i]
        du[i,IB] = b*dB[i] + ai*(1-b)*Bi*(Ai+2*Bi) +
                   (1-b)*(ai-1/ai)/Ji^2 + dBeta[i]/(ai*Ji)
    end
    return nothing
end
end # module EinsteinEquations
