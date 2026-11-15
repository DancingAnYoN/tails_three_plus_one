module ConstrainedHistory
using MMIEF: RunData
using GLMakie
using Printf
export plot_constraint_solved_history

"""Save a static four-panel observer/constraint history for a constraint-solved run."""
function plot_constraint_solved_history(data::RunData;
        outdir::AbstractString="results",show_window::Bool=false)
    GLMakie.activate!()
    mkpath(outdir)
    M=data.settings.M0
    T=(data.times .- data.shifts) ./ M
    t=data.times ./ M
    fig=Figure(size=(1500,1040),fontsize=17)
    Label(fig[0,1:2],"Constraint-solved Einstein–Klein–Gordon | MMIEF",
          fontsize=24,font=:bold)
    ax1=Axis(fig[1,1],xlabel="T / M₀",ylabel="φ",title="Full observer history")
    ax2=Axis(fig[1,2],xlabel="T / M₀",ylabel="log₁₀|φ|",title="Logarithmic amplitude")
    ax3=Axis(fig[2,1],xlabel="t / M₀",ylabel="Mass / M₀",title="Horizon and outer mass")
    ax4=Axis(fig[2,2],xlabel="t / M₀",ylabel="log₁₀|constraint|",title="Constraint diagnostics")
    palette=[:dodgerblue,:darkorange,:seagreen,:purple,:firebrick]
    for i in eachindex(data.settings.observers)
        J=data.settings.observers[i]
        y=data.signals[i]
        color=palette[mod1(i,length(palette))]
        name=@sprintf("J = %.0f M₀",J/M)
        lines!(ax1,T,y;color=color,linewidth=2,label=name)
        lines!(ax2,T,log10.(abs.(y).+1e-300);color=color,linewidth=2,label=name)
    end
    lines!(ax3,t,data.masses_h ./ M;label="Horizon",linewidth=2.4,color=:firebrick)
    lines!(ax3,t,data.masses_out ./ M;label="Outer Misner–Sharp mass",
           linestyle=:dash,linewidth=2.0,color=:slategray)
    lines!(ax4,t,log10.(abs.(data.ham_l2).+1e-25);
           label="Hamiltonian L₂",linewidth=2.0,color=:dodgerblue)
    lines!(ax4,t,log10.(abs.(data.scalar_constraint_l2).+1e-25);
           label="Φ−∂ᵣφ L₂",linewidth=2.0,color=:seagreen)
    lines!(ax4,t,log10.(abs.(data.horizon_error).+1e-25);
           label="Horizon residual",linewidth=2.0,color=:firebrick)
    axislegend(ax1,position=:rt)
    axislegend(ax2,position=:rt)
    axislegend(ax3,position=:rb)
    axislegend(ax4,position=:rt)
    Label(fig[3,1:2],"Constraint enforcement is not a proof of convergence; test multiple grids, boundary radii, and fit windows.",
          fontsize=14,tellwidth=false)
    path=joinpath(outdir,"constrained_history.png")
    save(path,fig)
    println("Saved full amplitude/constraint history: ",path)
    show_window && display(fig)
    return path
end
end # module ConstrainedHistory
