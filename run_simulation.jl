# VS Code / Julia REPL convenience entry point.
# Load once: include("run_simulation.jl")
# Then: main(["--smoke", "--vacuum"])
# Thin wrapper around scripts/run_simulation.jl (baseline solver).
include(joinpath(@__DIR__,"scripts","run_simulation.jl"))
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
