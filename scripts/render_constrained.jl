# Render a constraint-solved run: constrained_history.png and the standard
# dashboard, plus the pulse and dashboard movies with --movie.
# julia --project=viz scripts/render_constrained.jl results/constrained_quick
# julia --project=viz scripts/render_constrained.jl results/constrained_quick --movie
# Julia REPL: include("scripts/render_constrained.jl"); constrained_render_main(["results/constrained_quick","--show"])
# Visualization modules are included at top level to avoid world-age issues.
import MMIEF
include(joinpath(@__DIR__,"..","viz","Visualization.jl"))
import .MMIEFVisualization
include(joinpath(@__DIR__,"..","viz","ConstrainedHistory.jl"))
import .ConstrainedHistory

# args: run directory, then optional --movie and --show.
function constrained_render_main(args::Vector{String}=ARGS)
    isempty(args) && error("Usage: render_constrained.jl results/constrained_quick [--movie] [--show]")
    allowed=Set(["--movie","--show"])
    all(x->x in allowed,args[2:end]) || error("Unknown visualization option")
    outdir=abspath(args[1])
    data=MMIEF.load_run(outdir)
    show="--show" in args
    ConstrainedHistory.plot_constraint_solved_history(data;outdir=outdir,show_window=show)
    MMIEFVisualization.make_dashboard(data;outdir=outdir,movie="--movie" in args,
                                       show_window=false)
    if "--movie" in args
        MMIEFVisualization.make_pulse_movie(data;outdir=outdir,
                                             show_window=false)
    end
    return outdir
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    constrained_render_main(ARGS)
end
