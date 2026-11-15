# Constraint-solved (3+1) MMIEF Einstein–Klein–Gordon evolution

**Status:** new implementation / numerical research prototype. The existing
`MMIEF.run_evolution` is left intact; this extension adds
`MMIEF.run_constraint_solved`. Results are *not* convergence-validated.

## 1. Why this is more constrained than the previous formulation

The previous implementation evolves `a` and `B` using two Einstein equations,
then projects `a` onto the Hamiltonian constraint at intermediate RK4 stages.
This alternative **never advances `a` by a differential evolution equation**.
Instead, every evaluation of the method-of-lines right-hand side uses a newly
reconstructed spatial metric satisfying the radial Hamiltonian constraint.

The momentum constraint determines `A` from the derivative of `B`. Thus **both
independent Einstein constraints are enforced at every RK4 stage**, modulo
spatial discretization and the approximate outer scalar boundary condition.
The independent Einstein equation for `a` can be evaluated as a *residual*
without time-evolving `a`.

**Terminological caveat:** This is a *constraint-solved gravitational
formulation*, but is **not** an ADM *fully constrained formulation* in the
strict sense of solving for *all* gravitational fields from constraints and
gauge elliptic/radial equations without any gravitational evolution. The
remaining angular extrinsic curvature `B` still evolves. A purely fully
constrained system in this sense requires an additional compatible geometric
condition or a different characteristic formulation. Calling a Hamiltonian
projection alone fully constrained would be misleading.

## 2. Conventions, fields, and gauge

- Signature `(-,+,+,+)`, MTW curvature sign, `G=c=1`.
- Induced metric `gamma_ab = g_ab + n_a n_b`; ADM `K_ab = -1/2 L_n gamma_ab`.
- Metric: `ds² = -alpha² dt² + a²(dr + beta dt)² + J² dOmega²`.
- Shifted areal radius: `J(t,r)=r+f(t)`.
- Null slicing: `alpha = a(1-beta)`, so `v=t+r` is incoming null.
- Curvatures `A=K^r_r`, `B=K^theta_theta`.
- Scalar first-order variables `Phi = d_r phi`, `Pi = (a/alpha)(d_t phi-beta Phi)`.
- State `u[:, :] = (a, B, Phi, Pi, phi)` plus scalar `f`; column `a` is
  *constrained data* and has zero RK4 right-hand side.

From angular metric evolution and null slicing,

```
beta = (f_dot + a*J*B)/(1 + a*J*B)
alpha = a*(1 - beta)
```

The apparent horizon is at the fixed coordinate `r_h = 2*M0`; its areal
radius is `J_h = r_h + f`. Its condition `a_h*J_h*B_h=1` is enforced via the
Misner–Sharp mass boundary `m_h=J_h/2` with positive `B_h`.

The horizon-locking ODE is

```
f_dot = 4*pi*J_h^2/a_h^2 * (Phi_h + Pi_h)^2.
```

## 3. Hamiltonian and momentum constraints

The radial mass and radial curvature constraints are

```
m' = 4*pi*J^2 * [ (Phi^2+Pi^2)/(2*a^2) + (J*B/a)*Phi*Pi ]
a^(-2) = 1 - 2*m/J + (J*B)^2
A = B + J*B' - 4*pi*J/a * Phi*Pi
```

The integrator starts from `m(J_h)=J_h/2` and uses radial RK4 with the
repository's fourth-order midpoint interpolation. This reconstructs `a`
from the mass constraint throughout the slice, then reconstructs `A`
pointwise from the momentum constraint and differentiates other fields with
fourth-order finite differences.

The equation for `a_t` is **not** used to evolve `a`. The discretized
Hamiltonian and scalar reduction residuals are recorded in `diagnostics.csv`.
The momentum residual is algebraically zero relative to the chosen discrete
`B'`, although the reconstructed `A` is not written to a CSV column.

## 4. Evolved Einstein and Klein–Gordon equations

```
B_t = beta B' + a(1-beta) B(A+2B)
      + (1-beta)(a-1/a)/J^2 + beta'/(a J)

phi_t = beta Phi + (1-beta) Pi
Phi_t = [beta Phi+(1-beta)Pi]'
Pi_t = [J^2 (beta Pi+(1-beta)Phi)]'/J^2 - 2*f_dot*Pi/J
f_t = f_dot
```

The right hand side is evaluated at the *constrained* RK4 stage state.
The original stationary Schwarzschild discrete `B_t` residual is subtracted
as a well-balancing correction (the vacuum reference remains exactly static
at `f=0`). The incoming scalar characteristic at the outer boundary is
projected onto an approximate Sommerfeld condition iteratively with the
radial mass integration. This is an approximation at finite `rmax`.

## 5. Code reuse from the GitHub Actions-tested project

Direct reuse through `MMIEF` submodules:

| Repository module | Reused for |
|---|---|
| `Types.jl`, `Grid.jl` | Configuration and fixed radial grid |
| `FiniteDifferences.jl` | Fourth-order radial derivatives and KO filtering |
| `Integrators.jl`, `MethodOfLines.jl` | RK4 stages and CFL |
| `Constraints.jl` | Misner–Sharp radial mass projection, momentum reconstruction |
| `Geometry.jl` | Vacuum Schwarzschild and horizon locking |
| `ScalarField.jl` | Scalar flux RHS and Gaussian pulse |
| `InitialData.jl` | Coupled initial constraint data |
| `BoundaryConditions.jl` | Outer Sommerfeld projection |
| `DataOutput.jl`, `Observers.jl`, `Diagnostics.jl` | Standard CSV/TOML data, observers, constraints |
| `viz/Dashboard.jl`, `viz/PulseMovie.jl` | GLMakie dashboard and animation |

No existing numerical source module is deleted or overwritten except
`src/MMIEF.jl` gains a new `include`/export. The old `run_simulation.jl`
continues to select the previous solver.

## 6. Run the simulation in VS Code / Julia 1.12

From the root of this repository:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using MMIEF
include("scripts/run_constrained.jl")
constrained_main(["--smoke", "--vacuum"])
constrained_main(["--quick"])
# When the quick test passes:
constrained_main(["--standard"])
```

Or PowerShell:

```powershell
julia --project=. scripts/run_constrained.jl --smoke --vacuum
julia --project=. scripts/run_constrained.jl --quick
julia --project=. scripts/run_constrained.jl --standard
```

Output locations:

```
results/constrained_smoke_vacuum/
results/constrained_quick/
results/constrained_standard/
```

Each contains standard `diagnostics.csv`, `observers.csv`, `snapshots.csv`,
`final_state.csv`, `metadata.toml`, and `config.toml`. `metadata.toml`
identifies the constraint-solved evolution method.

## 7. GLMakie visualization

In the root of the repository (one-time setup):

```julia
using Pkg
Pkg.activate("viz")
Pkg.develop(path=abspath("."))
Pkg.instantiate()
```

From the Julia REPL:

```julia
include("scripts/render_constrained.jl")
constrained_render_main(["results/constrained_quick", "--show"])
constrained_render_main(["results/constrained_quick", "--movie"])
```

PowerShell:

```powershell
julia --project=viz scripts/render_constrained.jl results/constrained_quick --movie
```

The renderer produces

- `constrained_history.png`: entire signed/log-amplitude histories at fixed
  areal radii `J=10M0`, `30M0`, `50M0`, mass evolution, constraint residuals;
- `dashboard.png`: existing six-panel diagnostics, including exploratory QNM;
- with `--movie`: `dashboard.mp4` and `scalar_pulse_absorption.mp4`.

The visualization process reads saved CSVs and never reruns the PDE solver.
The visualization environment requires functional OpenGL and FFmpeg support.

## 8. Validation ladder / scientifically important limitations

1. **Vacuum:** `MMIEF.vacuum_constrained_rhs()` should be near machine zero;
   test `--smoke --vacuum` and longer vacuum runs.
2. **Constraint audit:** compare Hamiltonian L2, momentum reconstruction,
   horizon residual, and scalar reduction constraint across resolutions.
3. **Convergence:** run `dr=0.5,0.25,0.125` at matched physical domains and
   CFL; check fourth-order convergence in the smooth interior, not only
   single-run constraint smallness.
4. **Boundary:** repeat with larger `rmax` to determine the onset of numerical
   reflections. These are critical for late-time tails.
5. **Fixed background:** test scalar evolution against the analytic
   Schwarzschild potential; do not infer correct QNMs from horizon behavior.
6. **Mode extraction:** scalar `l=0` Schwarzschild theoretical mode
   `Mf*omega=0.110455-0.104896im`. The period (~56.9 Mf) is much larger
   than the e-folding time (~9.53 Mf); many ad hoc fit windows are unreliable.
7. **Tail:** expected generic `t^-3` fixed-radius scalar tail in 3+1 for
   weak linear scattering. Fit window, radius and domain size matter.

As of this package, no Julia runtime was available in the creation
container. **The new source has not been executed under Julia or GLMakie.**
Do not interpret generated astrophysical results until `Pkg.test()` and
multi-resolution numerical convergence tests pass on your machine.

## References

- R. Marsa and M. Choptuik, *Black-hole–scalar-field interactions in spherical
  symmetry*, Phys. Rev. D **54** (1996) 4929, arXiv:gr-qc/9607034.
- `EQUATIONS.md` and the recent 3+1 LaTeX MMIEF derivation notes.
