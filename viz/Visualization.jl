# GLMakie front end, loaded only in the separate viz/ environment so the core
# MMIEF package never depends on GLMakie. Scripts include this file at top level
# (not inside a function) to avoid Julia world-age errors.
module MMIEFVisualization
using MMIEF
using GLMakie
include("Dashboard.jl")
include("PulseMovie.jl")
using .Dashboard: make_dashboard
using .PulseMovie: make_pulse_movie
export make_dashboard, make_pulse_movie
end # module MMIEFVisualization
