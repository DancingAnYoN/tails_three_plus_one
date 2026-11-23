# MMIEF: spherically symmetric, massless Einstein–Klein–Gordon (EKG) evolution
# in minimally modified ingoing Eddington–Finkelstein (MMIEF) coordinates,
# after Marsa & Choptuik, Phys. Rev. D 54, 4929 (1996).
#
# Metric and gauge (see EQUATIONS.md):
#   ds² = -α² dt² + a² (dr + β dt)² + J² dΩ²,   J = r + f(t),   α = a(1-β).
#   r is the fixed numerical grid coordinate; J is the areal radius.
#   A = K^r_r, B = K^θ_θ with K_ij = -½ L_n γ_ij;  G = c = 1.
#   Scalar: Φ = ∂_r φ,  Π = (a/α)(∂_t φ - β Φ).
#
# Two solver paths share the physics modules below:
#   run_evolution          baseline ("partially constrained"): a and B are
#                          evolved; a is re-projected onto the Hamiltonian
#                          constraint at every RK4 stage.
#   run_constraint_solved  a is reconstructed from the Hamiltonian constraint
#                          at every stage and never evolved; B is evolved
#                          (see CONSTRAINED_FORMULATION.md).
# In both, A comes from the momentum constraint and f from horizon locking.
# Research prototype: not convergence-validated (see VALIDATION.md).
module MMIEF
# Include every submodule once, at package load time (Julia 1.12 world-age safe).
# Order matters: a module may only `using ..X` modules included before it.
include("Types.jl")              # Config, state column indices
include("Grid.jl")               # fixed radial grid r_i, config validation
include("FiniteDifferences.jl")  # 4th-order d/dr, Kreiss–Oliger filter
include("Geometry.jl")           # gauge (α, β), Schwarzschild data, horizon locking
include("ScalarField.jl")        # pulse profile, Klein–Gordon RHS
include("EinsteinEquations.jl")  # evolution equations for a and B
include("Constraints.jl")        # momentum (A) and Hamiltonian (a) constraint solves
include("InitialData.jl")        # constraint-satisfying initial slice
include("BoundaryConditions.jl") # approximate outgoing conditions at r_max
include("Integrators.jl")        # generic RK4 with projection hooks
include("MethodOfLines.jl")      # PDE-agnostic MoL step and CFL step size
include("Observers.jl")          # sampling φ at fixed areal radius J_obs
include("Diagnostics.jl")        # masses and constraint residuals
include("DataOutput.jl")         # RunData, versioned CSV/TOML output
include("Evolution.jl")          # baseline solver path
include("Analysis.jl")           # exploratory QNM fit and local power index
include("ConstraintSolvedEvolution.jl")  # constraint-solved solver path

using .Types: Config, SimulationState
using .Grid: radial_grid, RadialGrid
using .Geometry: schwarzschild_fields, mmief_outgoing_delay
using .InitialData: initial_data
using .Diagnostics: diagnostics
using .DataOutput: RunData, save_run, load_run, save_csv, save_waveforms, waveform_filename
using .Evolution: run_evolution, vacuum_test
using .Analysis: fit_qnm, local_power_index, local_loglog_slope, fit_tail_model, fit_qnm_tail
using .ConstraintSolvedEvolution: run_constraint_solved, vacuum_constrained_rhs

export Config, SimulationState, RadialGrid, radial_grid, schwarzschild_fields,
       initial_data, diagnostics, run_evolution, vacuum_test, RunData,
       save_run, load_run, save_csv, save_waveforms, waveform_filename,
       fit_qnm, local_power_index, local_loglog_slope, fit_tail_model,
       fit_qnm_tail, mmief_outgoing_delay,
       run_constraint_solved, vacuum_constrained_rhs
end # module MMIEF
