# Installing the constraint-solved solver in an existing GitHub clone

This ZIP is based on the previously supplied project source plus its agent
skills. A `src/ConstraintSolvedEvolution.jl` submodule has been added, and
`src/MMIEF.jl` includes/exports it. The existing evolution remains available.

## Full repository ZIP

Unpack the full archive into a new folder and run the commands in
`CONSTRAINED_FORMULATION.md`. It does not contain `.git` or simulation results.

## Overlay patch ZIP

If your own GitHub clone has newer source or skill edits, prefer the overlay.
Extract it into your repository root, **then make the following two edits** to
existing files (unless already done):

1. At the end of `src/MMIEF.jl`, alongside other `include` declarations,
   put `include("ConstraintSolvedEvolution.jl")` after `include("Analysis.jl")`.
   Add `using .ConstraintSolvedEvolution: run_constraint_solved, vacuum_constrained_rhs`
   alongside other `using .Module` imports, and export both names.
2. In `test/runtests.jl`, add
   `include("integration/test_constrained_evolution.jl")` after the other
   integration tests.

If desired, add this action step to `.github/workflows/ci.yml`:

```yaml
      - name: Constraint-solved vacuum smoke evolution
        shell: bash
        run: julia --project=. scripts/run_constrained.jl --smoke --vacuum
```

The source module reuses public submodules in MMIEF from the repository
snapshot. If those APIs have changed on your branch, adapt the imports.

## Run

From a Julia REPL at the repository root:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using MMIEF
include("scripts/run_constrained.jl")
constrained_main(["--smoke","--vacuum"])
constrained_main(["--standard"])
```

For plots, activate `viz`, develop the local `MMIEF` package and load:

```julia
Pkg.activate("viz")
Pkg.develop(path=abspath("."))
Pkg.instantiate()
include("scripts/render_constrained.jl")
constrained_render_main(["results/constrained_standard", "--movie"])
```

## Scientific caveat

**This implementation solves both independent constraints at every RK4
stage, but still evolves the angular curvature `B`.** It is a
constraint-solved (not strict elliptic fully constrained) approach,
consistent with the MMIEF gauge and existing source modules. See the
main formulation document for details.
