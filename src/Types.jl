# Run configuration and the column layout of the method-of-lines state.
module Types
export Config, IA, IB, IPHI, IPI, IFIELD, SimulationState, fieldview

# Run configuration. Lengths and times are in geometrized units (G = c = 1);
# the defaults assume M0 = 1, so they can also be read in units of M0.
Base.@kwdef struct Config
    M0::Float64 = 1.0              # initial black-hole mass; horizon held at r_h = 2*M0
    rmax::Float64 = 600.0          # outer boundary in the grid coordinate r
    dr::Float64 = 0.25             # radial spacing; (rmax - 2*M0)/dr must be an integer
    cfl::Float64 = 0.25            # Courant number C = dt/dr (reduced if speeds exceed 1)
    tmax::Float64 = 400.0          # final coordinate time t
    pulse_r0::Float64 = 20.0       # centre r0 of the Gaussian scalar pulse
    pulse_sigma::Float64 = 2.0     # width σ in exp[-(r - r0)²/σ²]
    pulse_amplitude::Float64 = 1e-4  # amplitude ε; 0 gives a vacuum Schwarzschild run
    cutoff_start::Float64 = 6.0    # pulse is exactly zero for r <= cutoff_start ...
    cutoff_end::Float64 = 9.0      # ... and unmodified for r >= cutoff_end
    ko_strength::Float64 = 0.06    # Kreiss–Oliger dissipation coefficient σ_KO
    output_dt::Float64 = 0.5       # cadence of diagnostics and observer samples
    snapshot_dt::Float64 = 2.0     # cadence of radial profile snapshots
    profile_stride::Int = 4        # snapshots keep every profile_stride-th grid point
    observers::Vector{Float64} = [10.0, 30.0, 50.0]  # fixed areal radii J_obs (absolute values)
    # Fail early when the horizon mass exceeds the initial outer mass by
    # far more than can be attributed to truncation error. A weak, isolated
    # ingoing pulse cannot increase the BH mass by many solar/geometric units.
    check_mass_budget::Bool = true
    mass_budget_fraction::Float64 = 0.01
end

# Column indices of the state matrix u[i, column], row i = grid point r_i:
# radial metric a, angular curvature B = K^θ_θ, Φ = ∂_r φ, Π, and φ itself.
const IA, IB, IPHI, IPI, IFIELD = 1, 2, 3, 4, 5

"""Full method-of-lines state: five radial fields and scalar horizon shift."""
mutable struct SimulationState
    fields::Matrix{Float64}   # N × 5 matrix, columns indexed by IA ... IFIELD
    f::Float64                # horizon shift f(t); areal radius J = r + f
end
"""Read-only named field-column accessor (without copying)."""
function fieldview(state::SimulationState, name::Symbol)
    index = Dict(:a=>IA, :B=>IB, :Phi=>IPHI, :Pi=>IPI, :phi=>IFIELD)
    haskey(index,name) || throw(ArgumentError("Unknown field $name"))
    return view(state.fields,:,index[name])
end
end # module Types
