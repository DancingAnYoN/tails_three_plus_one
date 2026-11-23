# Render the six-panel diagnostics dashboard (PNG, optionally MP4) from a
# saved run directory.
# Requires visualization environment; see README.
#   julia --project=viz scripts/render_dashboard.jl results/quick
# Reads the saved CSV/TOML output only; the solver is not rerun.
import MMIEF
include(joinpath(@__DIR__,"..","viz","Visualization.jl"))
import .MMIEFVisualization

# args: run directory, then optional --movie (also write dashboard.mp4) and
# --show (open a window).
function render_main(args::Vector{String}=ARGS)
    isempty(args) && error("Usage: render_dashboard.jl <run_directory> [--movie] [--show]")
    run_dir=abspath(args[1])
    allowed=Set(["--movie","--show"])
    all(x->x in allowed,args[2:end]) || error("Unknown render option")
    run=MMIEF.load_run(run_dir)
    return MMIEFVisualization.make_dashboard(run;outdir=run_dir,
                    movie=("--movie" in args),show_window=("--show" in args))
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    render_main(ARGS)
end
