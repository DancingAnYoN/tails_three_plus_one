# Plot the full stored waveform history at every fixed-areal-radius observer
# in observers.csv (raw, log10, absolute, running envelope).
#   julia --project=viz scripts/render_observer_history.jl results/standard [--show]
# Writes observer_history_J<J1>_J<J2>....png into the run directory, e.g.
# observer_history_J10_J30_J50.png for the default observers. Standalone:
# reads the CSV directly and does not load MMIEF.
using GLMakie
using Printf

# -------------------------
# Lightweight CSV utilities
# -------------------------
# Split one CSV line into trimmed String fields (no quoting support).
function _split_csv_line(line::String)
    return String.(strip.(split(chomp(line), ',')))
end

# Read a CSV into Dict(column name => column of strings). Short rows are padded
# with "" and long rows truncated to the header length.
function read_simple_csv(path::AbstractString)
    isfile(path) || error("CSV file not found: $path")
    lines = readlines(path)
    isempty(lines) && error("Empty CSV file: $path")
    header = _split_csv_line(lines[1])
    cols = Dict{String,Vector{String}}(h => String[] for h in header)
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = _split_csv_line(line)
        if length(vals) < length(header)
            append!(vals, fill("", length(header) - length(vals)))
        elseif length(vals) > length(header)
            vals = vals[1:length(header)]
        end
        for (h, v) in zip(header, vals)
            push!(cols[h], v)
        end
    end
    return cols
end

# Return the first column whose name matches one of `aliases`, compared
# case-insensitively and tried in the order given; `nothing` if none match.
function find_col(cols::Dict{String,Vector{String}}, aliases::Vector{String})
    lower_map = Dict(lowercase(k) => k for k in keys(cols))
    for a in aliases
        key = get(lower_map, lowercase(a), nothing)
        key === nothing || return key
    end
    return nothing
end

# Parse a field as Float64; unparsable or empty fields become NaN.
tofloat(s::AbstractString) = try parse(Float64, s) catch; NaN end

# (column name, parsed values) for the first matching alias, or (nothing, []).
function get_num_col(cols::Dict{String,Vector{String}}, aliases::Vector{String})
    key = find_col(cols, aliases)
    key === nothing && return nothing, Float64[]
    return key, [tofloat(v) for v in cols[key]]
end

# -------------------------
# Observer data loading
# -------------------------
# Return Dict(J_obs => (time, φ)) for every observer radius found in
# observers.csv (long format: one row per sample and observer).
function load_observers(run_dir::AbstractString)
    path = joinpath(run_dir, "observers.csv")
    cols = read_simple_csv(path)

    # Matching is case-insensitive and in alias order, so the observer-time
    # column T_asymptotic = t - f must come before "T", which would also match
    # the coordinate-time column "t". Plain t is only a fallback.
    _, T   = get_num_col(cols, ["T_asymptotic", "time_asymptotic", "T", "time", "t"])
    _, J   = get_num_col(cols, ["J_obs", "Jobs", "J", "J_over_M0"])
    _, phi = get_num_col(cols, ["phi", "varphi", "field", "scalar"])

    if isempty(T) || isempty(J) || isempty(phi)
        error("Could not identify T/J_obs/phi columns in observers.csv")
    end

    # One series per distinct observer radius in the file (the run's
    # Config.observers, e.g. chosen with --J in the run scripts).
    data = Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()
    for target in sort(unique(filter(isfinite, J)))
        idx = findall(k -> isfinite(J[k]) && abs(J[k] - target) < 1e-8, eachindex(J))
        data[target] = (T[idx], phi[idx])
    end
    return data
end

# -------------------------
# Plotting
# -------------------------
function make_observer_history(run_dir::AbstractString; show_window::Bool=false)
    data = load_observers(run_dir)
    isempty(data) && error("No observer data found in $(joinpath(run_dir, "observers.csv"))")
    radii = sort(collect(keys(data)))

    fig = Figure(size=(1600, 1200))
    ax1 = Axis(fig[1,1], xlabel="T / M0", ylabel="ϕ", title="Entire waveform history at fixed areal radii")
    ax2 = Axis(fig[1,2], xlabel="T / M0", ylabel="log10 |ϕ|", title="Amplitude history: log10 |ϕ(T)|")
    ax3 = Axis(fig[2,1], xlabel="T / M0", ylabel="|ϕ|", title="Absolute amplitude history")
    ax4 = Axis(fig[2,2], xlabel="T / M0", ylabel="cumulative max |ϕ|", title="Running envelope max{|ϕ|}")

    for J in radii
        t, phi = data[J]
        label = @sprintf("J = %g M0", J)
        lines!(ax1, t, phi, label=label)
        # 1e-300 avoids log10(0) at exact zeros.
        lines!(ax2, t, log10.(abs.(phi) .+ 1e-300), label=label)
        lines!(ax3, t, abs.(phi), label=label)
        env = accumulate(max, abs.(phi))
        lines!(ax4, t, env, label=label)
    end

    axislegend(ax1, position=:rb)
    axislegend(ax2, position=:rb)
    axislegend(ax3, position=:rb)
    axislegend(ax4, position=:rb)

    Label(fig[3,1:2],
          "This figure shows the full stored time history of the scalar amplitude at the fixed-areal-radius observers J = " *
          join([@sprintf("%gM₀", J) for J in radii], ", ") * ".\n" *
          "Top left: raw waveform. Top right: logarithmic amplitude. Bottom left: absolute amplitude. Bottom right: running envelope maximum.",
          tellwidth=false, justification=:left)

    outpath = joinpath(abspath(run_dir),
                       "observer_history_" * join([@sprintf("J%g", J) for J in radii], "_") * ".png")
    save(outpath, fig)
    println("Saved observer-history figure: $outpath")

    if show_window
        display(fig)
    end
    return outpath
end

function observer_history_main(args::Vector{String}=ARGS)
    isempty(args) && error("Usage: render_observer_history.jl <run_directory> [--show]")
    run_dir = abspath(args[1])
    allowed = Set(["--show"])
    all(x -> x in allowed, args[2:end]) || error("Unknown render option")
    return make_observer_history(run_dir; show_window=("--show" in args))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    observer_history_main(ARGS)
end
