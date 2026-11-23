# Spatial operators on the uniform radial grid. No physics in this module.
module FiniteDifferences
export deriv4!, ko6!
"""Fourth-order accurate 5-point derivative, including 4th-order one-sided endpoints."""
function deriv4!(d::AbstractVector{Float64}, y::AbstractVector{Float64}, h::Float64)
    n = length(y)
    n >= 7 || error("At least 7 grid points are required")
    z = 1.0 / (12h)
    @inbounds begin
        # Inner edge (horizon): fully one-sided at row 1, biased at row 2.
        # These rows also act as the inner "boundary condition": the horizon is
        # an outflow boundary, so no data are imposed there.
        d[1] = z*(-25y[1] + 48y[2] - 36y[3] + 16y[4] - 3y[5])
        d[2] = z*(-3y[1] - 10y[2] + 18y[3] - 6y[4] + y[5])
        # Interior: centred (y[i-2] - 8y[i-1] + 8y[i+1] - y[i+2]) / (12h).
        for i in 3:n-2
            d[i] = z*(y[i-2] - 8y[i-1] + 8y[i+1] - y[i+2])
        end
        # Outer edge: mirror images of the inner-edge formulas.
        d[n-1] = z*(3y[n] + 10y[n-1] - 18y[n-2] + 6y[n-3] - y[n-4])
        d[n] = z*(25y[n] - 48y[n-1] + 36y[n-2] - 16y[n-3] + 3y[n-4])
    end
    return d
end

"""Sixth-derivative KO filter: negative semidefinite for Fourier modes."""
function ko6!(du::Matrix{Float64}, u::Matrix{Float64}, h::Float64, strength::Float64)
    # Adds (σ/64h)·Δ⁶u to the time derivative du, where Δ⁶ is the 7-point sixth
    # difference (1, -6, 15, -20, 15, -6, 1). Its Fourier symbol is
    # -(4 sin²(kh/2))³ ≤ 0, so the term damps grid-scale modes; for smooth data
    # it is ≈ (σ/64) h⁵ ∂⁶u, below the fourth-order truncation error. It acts on
    # every column but skips the three points nearest each boundary.
    strength == 0 && return
    n = size(u, 1)
    weight = strength / (64h)
    @inbounds for col in 1:size(u, 2), i in 4:n-3
        du[i,col] += weight*(u[i-3,col] - 6u[i-2,col] + 15u[i-1,col]
                            - 20u[i,col] + 15u[i+1,col] - 6u[i+2,col] + u[i+3,col])
    end
end

end # module FiniteDifferences
