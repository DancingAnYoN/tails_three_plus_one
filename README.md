# MMIEF — modular Julia Einstein–Klein–Gordon solver

This is a refactor of the existing MMIEF partially constrained research
prototype, **not** a new verified high-order solver. Geometry and physics
follow Marsa–Choptuik; the time integrator is Method of Lines + RK4 with
Hamiltonian projection at intermediate stages and Schwarzschild discrete
residual subtraction. Ingoing null slicing: `alpha=a*(1-beta)`; shifted
areal radius: `J=r+f(t)`.

## 0. Quickest route: one run for a chosen (T, J)

```powershell
julia --project=. scripts/run_TJ.jl 600 30          # T = 600 M0, J = 30 M0
julia --project=. scripts/run_TJ.jl 600 10,30,50    # several observers
julia --project=. scripts/run_TJ.jl 600 30 --dr 0.5 --no-plot
```

REPL: `include("scripts/run_TJ.jl"); run_TJ(600, 30)`.

Fit the tail `phi = A t^-gamma exp(B/t + C/t^2)` of a finished run (linear
least squares in `ln|phi|`; also prints simpler models and the dependence on
the window start, and writes `tail_fit.csv`, `tail_fit_J<J>.csv`):

```powershell
julia --project=. scripts/fit_tail.jl results/T1000_J30_run --tmin 250
```

This runs the constraint-solved scheme (`a` and `A` from the constraints at
every RK4 stage, `B` evolved; see `CONSTRAINED_FORMULATION.md`) up to
coordinate time `t = T`. The outer boundary is placed automatically at
`rmax >= T + J + 10`, so no boundary signal can reach the observer. Results
go to `results/T600_J30_run/` (`results/T600_J10-30-50_run/`,
`results/T600_J30_dr0.5_run/` for the other examples), including
`phi_J30.csv` with `(t, varphi(t,J))` sampled every `0.5 M0` and the tail
plot `tail_loglog.png` (needs the `viz/` environment, section 2). Rerunning
the same parameters overwrites the run files in that folder; other files you
put there are kept.

## 1. Install and run the numerical solver (no GLMakie needed)

Open the folder `MMIEF_Modular` in VS Code, start a Julia REPL, and type:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using MMIEF
include("run_simulation.jl")
main(["--smoke", "--vacuum"])
main(["--quick"])
```

Or from PowerShell inside the project root:

```powershell
julia --project=. scripts/run_simulation.jl --smoke --vacuum
julia --project=. scripts/run_simulation.jl --quick
```

Available simulation sizes are `--smoke`, `--quick`, `--standard`, `--long`;
add `--vacuum` for Schwarzschild vacuum. Results go into
`results/smoke_vacuum/`, `results/quick/`, etc.

On top of a size you can choose the observers and the run length (units of
`M0`; both options also work with `scripts/run_constrained.jl`):

```powershell
julia --project=. scripts/run_simulation.jl --standard --J 30 --tmax 600
julia --project=. scripts/run_simulation.jl --long --J 10,30,50,100 --tmax=1200
```

- `--J J1[,J2,...]`: fixed areal radii `J_obs` at which `varphi` is sampled
  (default `10,30,50`); each must satisfy `2M0 < J_obs < rmax` of the preset.
- `--tmax T`: final **coordinate** time `t_max` (default: the preset value).
- With either option the results go to a folder named after the run, e.g.
  `results/standard_T600_J30/` or `results/constrained_long_T1200_J10-30-50-100/`,
  so the preset folder `results/standard/` is never overwritten.
  Observer time `T=t-f` differs from `t` by `f=2(M_h-M0)`, which is tiny for
  the default weak pulse.

Each observer gets a two-column file `phi_J<J>.csv` (header `t,phi`) holding
`(t, varphi(t,J_obs))` at every output sample (cadence `output_dt=0.5`). The
script warns when `t_max > rmax - J_obs`: after that time the outer boundary
can influence the observer, so use a larger preset for late-time tails.
In the REPL: `main(["--standard", "--J", "30", "--tmax", "600"])`.

Add `--plot` to also write `tail_loglog.png` into the run directory: log–log
`|varphi|` against `T=t-f` with a `T^-3` guide, and the local power index
`p(T)=-d ln|varphi|/d ln T`. It runs `scripts/render_tail.jl` in a separate
process with the `viz/` environment (set up once, section 2); if that fails,
the run data are still saved. For an existing run directory:

```powershell
julia --project=viz scripts/render_tail.jl results/standard
```

**Do not** use `--no-movie` or `--no-gl` with the numerical runner: visualization
is now a separate process. This is intentional.

## 2. Install the separate GLMakie visualization environment

From the root folder in the Julia REPL:

```julia
using Pkg
Pkg.activate("viz")
Pkg.develop(path=abspath("."))
Pkg.instantiate()
```

If you use a fresh Julia REPL, the terminal equivalents are:

```powershell
julia --project=viz scripts/render_dashboard.jl results/quick
julia --project=viz scripts/render_movie.jl results/quick
```

These load **saved numerical results**; they do not rerun the PDE simulation.
Outputs include `dashboard.png` and `scalar_pulse_absorption.mp4`.
Add `--movie` to `render_dashboard.jl` to render the six-panel dashboard MP4.

When manually switching between the core and viz environments, activate the
right one with `Pkg.activate(".")` or `Pkg.activate("viz")`.

## 3. Tests

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

Tests are grouped as `test/unit`, `test/integration`, and `test/regression`.
Core tests do not load GLMakie. For a separate OpenGL-dependent smoke test:

```powershell
julia --project=viz test/visual_smoke.jl
```

## 4. Data format

Every run directory contains:

- `config.toml`: parameters used in the run.
- `metadata.toml`: schema v1, coordinate and unit conventions.
- `diagnostics.csv`: horizon mass, outer mass and constraint norms.
- `observers.csv`: long-format `(t,T_asymptotic,J_obs,phi)` time series.
- `phi_J<J>.csv`: one per observer, `(t,phi)` at that `J_obs` (same samples
  as `observers.csv`; not read by `load_run`).
- `snapshots.csv`: long-format `(t,f,r,J,Jphi)` movie input data.
- `final_state.csv`: complete five-column evolved state at final time.

You can read it from any Julia process with `MMIEF.load_run("results/quick")`.
All CSV values are plain text and intended to be portable. Large `--long`
runs may create sizeable snapshot files; adjust `snapshot_dt` and
`profile_stride` if disk use becomes too high. **There is not yet a restart
checkpoint format**: final_state.csv alone cannot resume integration.

## 5. Source layout

- `src/Types.jl`, `Grid.jl`: configuration, field indexes, grid.
- `src/FiniteDifferences.jl`, `Integrators.jl`, `MethodOfLines.jl`: independent MoL methods.
- `src/Geometry.jl`, `ScalarField.jl`, `EinsteinEquations.jl`:
  continuum formulas in semi-discrete RHS form.
- `src/Constraints.jl`, `InitialData.jl`, `BoundaryConditions.jl`:
  constraint solve/projection, Gaussian, excision & outgoing boundaries.
- `src/Evolution.jl`: coupled RHS and run driver.
- `src/Diagnostics.jl`, `Observers.jl`, `DataOutput.jl`:
  monitoring, fixed-areal-radius sampling, versioned CSV/TOML.
- `src/Analysis.jl`: exploratory QNM and tail estimators.
- `viz/Dashboard.jl`, `PulseMovie.jl`: GLMakie (optional environment).
- `scripts/`: standalone runner, renderer and analysis commands.
- `legacy/MMIEFCore.jl`: frozen monolithic solver used for equivalence tests.

## Validation status

Static source refactoring and independent Python checks were possible in the
creation environment; Julia and GLMakie were not installed there. The included
Julia tests **have not been executed here**. Run `Pkg.test()` and the optional
GLMakie smoke test locally before treating simulation output as reliable.
The original research prototype is also not convergence-validated for late
time Price tails or quasinormal modes. See `VALIDATION.md`.

## Constraint-solved MMIEF alternative (new)

A new branch of the numerical method enforces the Hamiltonian radial mass
constraint and the momentum constraint at **every RK4 stage**, reconstructing
`a` and `A` rather than time-evolving `a`. The remaining Einstein variable
`B` continues to evolve. See [CONSTRAINED_FORMULATION.md](CONSTRAINED_FORMULATION.md)
for the precise distinction between this formulation and a strictly fully
constrained formulation, equations, tests, and known limitations.

```powershell
julia --project=. scripts/run_constrained.jl --smoke --vacuum
julia --project=. scripts/run_constrained.jl --quick
julia --project=. scripts/run_constrained.jl --standard
julia --project=viz scripts/render_constrained.jl results/constrained_standard --movie
```

This addition leaves the original `run_simulation.jl` implementation intact.
