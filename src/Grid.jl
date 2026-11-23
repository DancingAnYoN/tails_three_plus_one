# Fixed, uniform radial grid r_i = 2*M0 + (i-1)*dr, i = 1..N.
# The first point is the apparent horizon r_h = 2*M0 (inner excision boundary).
# The grid never moves; the physical areal radius J_i(t) = r_i + f(t) shifts
# rigidly outward as the black hole grows.
module Grid
using ..Types: Config
export RadialGrid, check_config, radial_grid

struct RadialGrid
    r::Vector{Float64}
    dr::Float64
end
"""Validate a `Config` before a grid is built; throws on the first violated assumption."""
function check_config(p::Config)
    rh = 2*p.M0
    p.M0 > 0 || error("M0 must be positive")
    # Enough points for the 5-point stencils and the 7-point KO filter.
    p.dr > 0 && p.rmax > rh + 12*p.dr || error("Invalid radial grid")
    nsteps = (p.rmax-rh)/p.dr
    abs(nsteps-round(nsteps)) < 1e-8 ||
        error("(rmax-2M0)/dr must be an integer to place the outer boundary exactly")
    # Conservative bound for RK4 with these stencils and KO dissipation.
    0 < p.cfl <= 0.5 || error("Choose 0 < cfl <= 0.5 for this implementation")
    p.tmax > 0 && p.output_dt > 0 && p.snapshot_dt > 0 || error("Invalid time settings")
    p.pulse_sigma > 0 && p.pulse_amplitude >= 0 || error("Invalid pulse")
    # The cutoff keeps the scalar field exactly zero at the horizon initially.
    rh < p.cutoff_start < p.cutoff_end < p.pulse_r0 || error("Cutoff must lie between horizon and pulse")
    # Initial check only: observers can still leave the grid later as f grows,
    # which the evolution loops detect.
    all(J -> rh < J < p.rmax, p.observers) || error("Observers must lie in domain")
    0 <= p.ko_strength <= 0.15 || error("Use a small KO strength")
    p.profile_stride >= 1 || error("profile_stride must be positive")
    p.mass_budget_fraction >= 0 || error("mass_budget_fraction must be nonnegative")
end

"""Return the grid coordinates r_1 = 2*M0, ..., r_N = rmax after validating `p`."""
function radial_grid(p::Config)
    check_config(p)
    # check_config guarantees (rmax - 2M0)/dr is an integer up to 1e-8; round,
    # not floor, so rounding error (e.g. dr = 0.35) cannot drop the last point.
    N = round(Int, (p.rmax - 2*p.M0) / p.dr) + 1
    return [2*p.M0 + (i - 1) * p.dr for i in 1:N]
end

RadialGrid(p::Config) = RadialGrid(radial_grid(p), p.dr)
end # module Grid
