# Print a quick single-mode ringdown fit for every observer of a saved run.
#   julia --project=. scripts/analyze_waveforms.jl results/standard
# Exploratory only; scripts/check_qnm_vs_theory.jl scans several fit windows.
using MMIEF
using Printf

function analyze_main(args::Vector{String}=ARGS)
    length(args)==1 || error("Usage: analyze_waveforms.jl <run_directory>")
    data=load_run(abspath(args[1]))
    # Observer time T = t - f, and the final horizon mass as M_f.
    times=data.times.-data.shifts
    M=data.masses_h[end]
    println("Exploratory analysis only — validate convergence and fit windows.")
    for (j,J) in enumerate(data.settings.observers)
        phi=data.signals[j]
        @printf("J/M0 = %.2f; max |phi|=%.5e\n",J/data.settings.M0,maximum(abs,phi))
        # Fixed window T ∈ [45, 110] M0 for every observer (not delay-shifted).
        fit=fit_qnm(times,phi;tlo=45*data.settings.M0,
                    thi=min(110*data.settings.M0,times[end]),mass=M)
        if fit===nothing
            println("    QNM fit unavailable (short run).")
        else
            @printf("    exploratory M omega_R=%.6f, M gamma=%.6f, residual=%.3e\n",
                    fit.omegaM,fit.gammaM,fit.rel_error)
        end
    end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    analyze_main(ARGS)
end
