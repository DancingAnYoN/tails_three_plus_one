# Import this version into your already-cloned GitHub repository

This package was prepared from your uploaded archive `late-time_tails_in_einstein_klein_gordon_system.zip`.
The numerical source files are unchanged. The embedded `.git/` directory has been deliberately removed
so its old history/remotes cannot overwrite your new repository's Git configuration.

## 1. Copy the files into your empty clone (Windows PowerShell)

Suppose you already ran `git clone https://github.com/YOUR-USERNAME/YOUR-REPO.git`.
Extract this ZIP, then copy **the contents of this folder** into the clone.
Do not copy a second top-level folder into the clone, and do not use `git init` again.

For example, from PowerShell (edit both paths):

```powershell
$source = 'D:\Downloads\late-time_tails_in_einstein_klein_gordon_system'
$clone  = 'D:\PROJECT_PHYSICS\YOUR-REPO'
Get-ChildItem -LiteralPath $source -Force | Copy-Item -Destination $clone -Recurse -Force
Set-Location -LiteralPath $clone
```

The `-Force` ensures that hidden files such as `.gitignore` and `.github/` are copied.
Your cloned folder should retain its original `.git/` directory (which is not in this archive).

## 2. Verify and publish

```powershell
git status
git add .
git status --short
git commit -m "Import uploaded MMIEF late-time tails project"
git branch -M main
git push -u origin main
```

If `main` already has commits from GitHub, do not force-push. Use a normal `git pull` /
reconciliation workflow before pushing. Your existing `.git` history and `origin` remote remain intact.

The provided `.gitignore` ignores `results/` output and `Manifest.toml` files.
If using the working-copy ZIP, simulation results remain available locally but are not staged for upload.
You can confirm with `git status --short` before committing.

## 3. Run numerical and visualization checks

From the repository root:

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=. scripts/run_simulation.jl --smoke --vacuum
julia --project=. scripts/run_simulation.jl --standard
```

GLMakie visualization environment: open the Julia REPL from the repo root, then run:

```julia
using Pkg
Pkg.activate("viz")
Pkg.develop(path=abspath("."))
Pkg.instantiate()
```

From PowerShell, after generating or copying local `results/standard`:

```powershell
julia --project=. scripts/check_qnm_vs_theory.jl results/standard
julia --project=viz scripts/render_qnm_visualization.jl results/standard
julia --project=viz scripts/render_observer_history.jl results/standard
```

The QNM analysis is exploratory. Check convergence and fitting-window sensitivity before drawing scientific conclusions.

## Files preserved and excluded

- Preserved: all `src/`, `scripts/`, `viz/`, `test/`, `legacy/`, documentation, license, and GitHub workflow files.
- Excluded: the uploaded `.git/` directory, because the clone already contains correct Git metadata.
- Source-only archive: also excludes local `results/` outputs and package `Manifest.toml` files.
- Working-copy archive: retains local `results/` and `Manifest.toml` files; both are ignored by Git.
