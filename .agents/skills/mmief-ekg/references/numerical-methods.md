# MMIEF numerical implementation and commands

## The numerical state

Store `u[i,:]=(a,B,Phi,Pi,varphi)` at fixed uniform grid points
`r_i=r_h+i Delta r`; evolve scalar gauge variable `f(t)` separately.
At each RHS evaluation compute `J=r+f`, `fdot` from horizon flux,
`A` via radial momentum, `beta/alpha` via gauge, and scalar flux arrays.

Spatial derivative (fourth-order interior):

$$ (D_r u)_i=\frac{u_{i-2}-8u_{i-1}+8u_{i+1}-u_{i+2}}
{12\Delta r}. $$

Classical RK4 for state `Y=(u,f)`:

$$k_1=F(Y_n),\ \ k_2=F(Y_n+\Delta t k_1/2),\ \ 
 k_3=F(Y_n+\Delta t k_2/2),\ \ k_4=F(Y_n+\Delta t k_3),$$

$$Y_{n+1}=Y_n+\frac{\Delta t}{6}(k_1+2k_2+2k_3+k_4).$$

Here stage boundary conditions and Hamiltonian projection are explicit
callbacks; the generic RK4 integrator should remain physics-agnostic.
The physical step uses Courant factor `C=Delta t/Delta r` and reduces time
step if the observed characteristic-speed bound requires it.

The mass-constraint projection integrates `m'` radially from
`m(r_h)=J_h/2`, then reconstructs `a`. This intentionally modifies
stage states and requires convergence testing. Vacuum well balancing subtracts
a precomputed discrete Schwarzschild gravitational RHS to limit drift.
An optional KO6 term damps high-frequency grid modes.

## Numerical modes (from scripts/run_simulation.jl)

| Mode | rmax/M0 | Delta r/M0 | tmax/M0 | output Delta t/M0 |
|---|---:|---:|---:|---:|
| `--smoke` | 70 | 0.5 | 2 | 0.5 |
| `--quick` | 160 | 0.5 | 65 | 0.5 |
| `--standard` | 700 | 0.25 | 350 | 0.5 |
| `--long` | 1700 | 0.25 | 800 | 0.5 |

Use `--vacuum` for zero scalar amplitude (vacuum baseline).
`--J J1[,J2,...]` replaces the observer radii and `--tmax T` the final
coordinate time (units of M0; both run scripts, via `scripts/run_options.jl`);
such runs are saved to `results/<mode>[_vacuum]_T<T>_J<J>/`, never to the
preset folder.
The scripts warn when `tmax > rmax - J_obs` (outer boundary in the observer's
causal past). `scripts/run_TJ.jl T J [--dr DR] [--no-plot]` runs the
constraint-solved scheme with `rmax >= T + J + 10`, Courant number 0.125, and
writes `results/T<T>_J<J>[_dr<DR>]_run/`.

## Commands (shell, project root)

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=. scripts/run_simulation.jl --smoke --vacuum
julia --project=. scripts/run_simulation.jl --standard
julia --project=. scripts/check_qnm_vs_theory.jl results/standard
```

One-time visualization setup from Julia REPL in project root:

```julia
using Pkg
Pkg.activate("viz")
Pkg.develop(path=abspath("."))
Pkg.instantiate()
```

Graphics and movie (shell):

```powershell
julia --project=viz scripts/render_dashboard.jl results/standard
julia --project=viz scripts/render_movie.jl results/standard
julia --project=viz scripts/render_observer_history.jl results/standard
julia --project=viz scripts/render_qnm_visualization.jl results/standard
julia --project=viz scripts/render_tail.jl results/standard
```

`render_tail.jl` (also run by `--plot` in both run scripts) writes
`tail_loglog.png`: log–log |φ| vs T = t − f with a T⁻³ guide, and p(T) from
`local_loglog_slope` (direct local fit, only where φ keeps one sign). Do not
use the peak-envelope `local_power_index` for monotone tails without
correcting its ≈ T/(T−h) upward bias.

`scripts/fit_tail.jl <run> [--tmin 150]` fits `A t^-γ exp(b1/t + b2/t²)`
(`MMIEF.fit_tail_model`, linear least squares in ln|φ|) on
`[tmin, min(t_end, rmax − J)]`, with lower-order models and a window-start
scan; quote γ with its window spread, and check resolution separately.

`results/standard/` must exist **before** analysis/visualization;
GitHub source archives omit generated simulation data.

REPL commands (if needed):

```julia
using Pkg
Pkg.activate(".")
using MMIEF
include("run_simulation.jl")
main(["--standard"])

Pkg.activate("viz")
include("scripts/render_movie.jl")
movie_main(["results/standard"])
```

Use a fresh REPL if old global function definitions conflict.

## Output and reproducibility

- `config.toml`, `metadata.toml`: parameters, units, coordinate conventions.
- `diagnostics.csv`: time, horizon and outer mass, constraint residuals.
- `observers.csv`: long-format `t,T_asymptotic,J_obs,phi`.
- `phi_J<J>.csv`: per observer, `t,phi` (copy of the observers.csv samples).
- `snapshots.csv`: `t,f,r,J,Jphi`.
- `final_state.csv`: final `u` (not a restart checkpoint).

For fixed areal radius `J_obs` use `r_obs=J_obs-f(t)` and interpolate;
never sample at a stationary `r=J_obs`. Report both `M0` and final `Mf`.

## Pitfalls from prior runs

- Julia parses `2f2` as a possible Float32 literal, not `2*f2`.
- An imported package `using MMIEF` does not define script-local `main`.
  `include("run_simulation.jl")` is needed before `main(...)`.
- The QNM script defines `qnm_check` only after `include`ing it.
- The dashboard script entry is `render_main`, not `dashboard_main`.
- `get_num_col(::Dict{String,...},...)` does not accept
  `Dict{SubString{String},...}`. Normalize headers to `String`.
- GLMakie belongs under `viz`, not the solver environment.
- User-level simulations can fail if the moving observer leaves the grid;
  do not silently clamp its coordinate to the nearest grid point.
