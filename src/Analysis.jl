# Exploratory waveform analysis: single-mode ringdown fit and local power
# index of late-time tails. Results depend on fit windows and must be checked
# across windows, observers and resolutions before physical interpretation.
module Analysis
export fit_qnm, local_power_index, local_loglog_slope, fit_tail_model, fit_qnm_tail
"""Exploratory damped-sinusoid fit by a (gamma,omega) grid search.
Not a robust multi-overtone or Prony/QNM estimator: inspect window dependence.

Model on the window [tlo, thi], with τ = t - tlo and no constant offset:

    y ≈ e^{-γτ}(C cos ωτ + S sin ωτ).

For each (γM, ωM) on a grid the amplitudes C, S are solved exactly by linear
least squares; the pair with the smallest residual wins. A 30×30 coarse grid is
followed by a 23×23 refinement around the best point. Pass `t` as the observer
time T = t - f and `mass` as the final mass, so γM and ωM are dimensionless.
"""
function fit_qnm(t::Vector{Float64}, y::Vector{Float64};
                 tlo::Float64=45.0, thi::Float64=100.0,
                 mass::Float64=1.0)
    inds=findall(i -> tlo<=t[i]<=thi, eachindex(t))
    length(inds)>=25 || return nothing
    tt=[t[i]-tlo for i in inds]
    yy=y[inds]
    best=(Inf, 0.0,0.0,0.0,0.0)
    # Dimensionless search in M omega; similar scalar l=0 Schwarzschild domain.
    for (gmin,gmax,wmin,wmax,nsteps) in ((0.045,0.16,0.055,0.18,30),
                                          (NaN,NaN,NaN,NaN,23))
        if isnan(gmin)
            _,g0,w0,_,_=best
            dg=0.115/(30-1)*2
            dw=0.125/(30-1)*2
            gmin=max(0.02,g0-2dg); gmax=g0+2dg
            wmin=max(0.02,w0-2dw); wmax=w0+2dw
        end
        for gammaM in range(gmin,gmax,length=nsteps), omegaM in range(wmin,wmax,length=nsteps)
            gamma=gammaM/mass; omega=omegaM/mass
            # Accumulate the 2×2 normal equations for basis functions
            # C(t) = e^{-γt}cos ωt and S(t) = e^{-γt}sin ωt.
            scc=0.0; sss=0.0; scs=0.0; syc=0.0; sys=0.0; syy=0.0
            for k in eachindex(tt)
                decay=exp(-gamma*tt[k])
                C=decay*cos(omega*tt[k]); S=decay*sin(omega*tt[k]); yk=yy[k]
                scc+=C*C; sss+=S*S; scs+=C*S
                syc+=yk*C; sys+=yk*S; syy+=yk*yk
            end
            det=scc*sss-scs*scs
            det<=1e-25 && continue
            ampC=(syc*sss-sys*scs)/det
            ampS=(sys*scc-syc*scs)/det
            # Least-squares residual ‖y - fit‖² = ‖y‖² - (ampC·⟨y,C⟩ + ampS·⟨y,S⟩).
            residual=max(0.0,syy-(ampC*syc+ampS*sys))
            if residual<best[1]
                best=(residual,gammaM,omegaM,ampC,ampS)
            end
        end
    end
    isfinite(best[1]) || return nothing
    res,gM,wM,C,S=best
    return (gammaM=gM, omegaM=wM, C=C, S=S, tlo=tlo, thi=thi,
            rel_error=sqrt(res/max(sum(abs2,yy),eps(Float64))), mass=mass)
end

"""Estimate local -dlog(envelope)/dlog(t) with a running peak envelope.
Window parameters are in coordinate time units. Avoid zero crossings by using
maximum |phi| over ±peak_halfwidth, then regress over ±fit_halfwidth.

Returns a vector of the same length as `t`, NaN where no estimate is made
(before `tstart`, below the amplitude floor, or with too few points). For a
generic l = 0 Price tail at fixed radius, φ ~ T⁻³, the index tends to 3.

Bias: for a monotonically decaying |y| the running maximum over ±h is
|y(t - h)|, so a pure power law t^-p gives about p t/(t - h) instead of p
(+4 % at t = 300, h = 12). Use `local_loglog_slope` for a non-oscillating tail.
"""
function local_power_index(t::Vector{Float64},y::Vector{Float64};
                           peak_halfwidth::Float64=15.0,
                           fit_halfwidth::Float64=25.0,
                           tstart::Float64=100.0)
    n=length(t)
    env=zeros(Float64,n)
    for i in 1:n
        l=searchsortedfirst(t,t[i]-peak_halfwidth)
        h=searchsortedlast(t,t[i]+peak_halfwidth)
        vmax=0.0
        @inbounds for j in l:h
            vmax=max(vmax,abs(y[j]))
        end
        env[i]=vmax
    end
    out=fill(NaN,n)
    flooramp=max(maximum(env)*1e-12,1e-29)
    for i in 1:n
        if t[i] < tstart || env[i] < flooramp
            continue
        end
        lo=searchsortedfirst(t,max(tstart,t[i]-fit_halfwidth))
        hi=searchsortedlast(t,t[i]+fit_halfwidth)
        hi-lo<12 && continue
        # Least-squares slope of log(env) against log(t) over the window.
        sx=0.0; sy=0.0; sxx=0.0; sxy=0.0; count=0
        for j in lo:hi
            if t[j]>0 && env[j]>=flooramp
                x=log(t[j]); z=log(env[j])
                sx+=x; sy+=z; sxx+=x*x; sxy+=x*z; count+=1
            end
        end
        den=count*sxx-sx*sx
        if count>12 && den>1e-12
            out[i]=-(count*sxy-sx*sy)/den
        end
    end
    return out
end

"""
    local_loglog_slope(t, y; halfwidth=25.0, tstart=0.0) -> Vector{Float64}

Local power index p(t) = -d ln|y| / d ln t, from a least-squares fit of ln|y|
against ln t over [t - halfwidth, t + halfwidth] (cut to [tstart, t_end]). No
envelope is used, so the estimate is exact for y = c t^-p. It is therefore only
defined where y keeps one sign over the whole window, i.e. on a non-oscillating
late-time tail, not during ringdown. NaN elsewhere, or with fewer than 5 points.
"""
function local_loglog_slope(t::Vector{Float64}, y::Vector{Float64};
                            halfwidth::Float64=25.0, tstart::Float64=0.0)
    n=length(t)
    length(y)==n || throw(DimensionMismatch("t and y must have equal length"))
    out=fill(NaN,n)
    for i in 1:n
        t[i] < tstart && continue
        lo=searchsortedfirst(t,max(tstart,t[i]-halfwidth))
        hi=searchsortedlast(t,t[i]+halfwidth)
        hi-lo < 4 && continue
        t[lo] > 0 || continue
        s=sign(y[lo])
        s != 0 || continue
        all(j -> sign(y[j])==s, lo:hi) || continue
        # Two-pass (centred) regression: ln t varies little across a late
        # window, so the one-pass form m*Σx² - (Σx)² would lose most digits.
        m=hi-lo+1
        xbar=0.0; zbar=0.0
        @inbounds for j in lo:hi
            xbar+=log(t[j]); zbar+=log(abs(y[j]))
        end
        xbar/=m; zbar/=m
        sxx=0.0; sxz=0.0
        @inbounds for j in lo:hi
            dx=log(t[j])-xbar
            sxx+=dx*dx; sxz+=dx*(log(abs(y[j]))-zbar)
        end
        sxx > 0 || continue
        out[i]=-sxz/sxx
    end
    return out
end

# Least-squares solution of min ‖X c - z‖ for a tall matrix with a few
# columns, by modified Gram–Schmidt QR with one re-orthogonalisation pass
# (keeps the conditioning of X, unlike the normal equations). Avoids a
# LinearAlgebra dependency in the core package.
function lstsq_small(X::Matrix{Float64}, z::Vector{Float64})
    m, n = size(X)
    m >= n || throw(ArgumentError("Need at least as many rows as columns"))
    Q = copy(X)
    R = zeros(n, n)
    for j in 1:n
        for _ in 1:2, k in 1:j-1
            r = sum(view(Q,:,k) .* view(Q,:,j))
            R[k,j] += r
            Q[:,j] .-= r .* view(Q,:,k)
        end
        R[j,j] = sqrt(sum(abs2, view(Q,:,j)))
        R[j,j] > 0 || error("Rank-deficient fit design matrix")
        Q[:,j] ./= R[j,j]
    end
    qz = [sum(view(Q,:,j) .* z) for j in 1:n]
    c = zeros(n)
    for j in n:-1:1
        s = qz[j]
        for k in j+1:n
            s -= R[j,k]*c[k]
        end
        c[j] = s/R[j,j]
    end
    return c
end

"""
    fit_tail_model(t, y; tmin, tmax=Inf, order=2) -> NamedTuple

Fit the late-time tail model

    y(t) = A t^(-γ) exp(b1/t + b2/t²)

to the samples with tmin <= t <= tmax (`order` = 1 drops b2, `order` = 0 is
the pure power law A t^-γ). In ln|y| the model is linear in (ln A, γ, b1, b2),
so this is an exact linear least-squares fit of ln|y|: every sample has the
same weight in *relative* error across the decades of decay, and no starting
values or iterations are needed. y must keep one sign on the window; A carries
that sign. b1/t + b2/t² absorbs the leading corrections to the asymptotic
power law (a time-origin offset c gives b1 ≈ γc, b2 ≈ γc²/2).

Returns (A, gamma, b1, b2, rms, maxabs, n, tmin, tmax); rms and maxabs are the
residuals of ln|y| (≈ relative residuals of y). Exploratory: the residuals are
systematic, not noise, so judge γ by its dependence on the window, observer
and resolution rather than by a formal error bar.
"""
function fit_tail_model(t::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                        tmin::Real, tmax::Real=Inf, order::Int=2)
    0 <= order <= 2 || throw(ArgumentError("order must be 0, 1 or 2"))
    length(t)==length(y) || throw(DimensionMismatch("t and y must have equal length"))
    idx = findall(i -> tmin <= t[i] <= tmax, eachindex(t))
    length(idx) >= order + 4 || error("Too few samples in the fit window")
    tw = Float64.(t[idx]); yw = Float64.(y[idx])
    tw[1] > 0 || error("The fit window must have t > 0")
    s = sign(yw[1])
    (s != 0 && all(v -> sign(v) == s, yw)) ||
        error("y changes sign or vanishes in the fit window; start the window later")
    # Columns scaled to O(1) with tref = first window time:
    # 1, ln(t/tref), tref/t, (tref/t)².
    tref = tw[1]
    X = hcat(ones(length(tw)), log.(tw ./ tref), [(tref ./ tw).^k for k in 1:order]...)
    z = log.(abs.(yw))
    c = lstsq_small(X, z)
    res = copy(z)                      # z - X c, column by column
    for j in eachindex(c)
        res .-= c[j] .* view(X,:,j)
    end
    gamma = -c[2]
    return (A = s*exp(c[1] + gamma*log(tref)), gamma = gamma,
            b1 = order >= 1 ? c[3]*tref : 0.0,
            b2 = order >= 2 ? c[4]*tref^2 : 0.0,
            rms = sqrt(sum(abs2, res)/length(res)), maxabs = maximum(abs, res),
            n = length(tw), tmin = tw[1], tmax = tw[end])
end

"""
    fit_qnm_tail(t, y; tlo, thi, tail=nothing, mass=1.0) -> NamedTuple or nothing

Single damped sinusoid plus an optional tail template on [tlo, thi]:

    y ≈ e^{-Γτ}(a cos ωτ + b sin ωτ) + κ tail,    τ = t - tlo,

where `tail` is a vector of template values at the samples `t` (or `nothing`
for no tail term). For each (Γ, ω) the amplitudes (a, b, κ) follow from linear
least squares (variable projection), and (MΓ, Mω) is found on a grid over
[0.03, 0.25]² (89 × 89), refined twice (41 × 41 over ±2 steps; final step
≈ 2.5e-5). Pass `t` as observer time T and `mass` as the final mass, so that
`omegaM` = M ω and `gammaM` = M Γ; to compare with M_f ω = 0.110455 - 0.104896 i
(scalar l = 0, exp(-iωT)).

Returns (omegaM, gammaM, a, b, kappa, rel_rms, tlo, thi, n, mass), with
kappa = NaN without a tail and rel_rms = ‖residual‖/‖y‖ on the window, or
`nothing` with fewer than 12 samples. Exploratory: compare windows, observers
and resolutions before identifying a QNM.
"""
function fit_qnm_tail(t::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                      tlo::Real, thi::Real, tail=nothing, mass::Real=1.0)
    length(t)==length(y) || throw(DimensionMismatch("t and y must have equal length"))
    tail === nothing || length(tail)==length(t) ||
        throw(DimensionMismatch("tail template must have one value per sample"))
    idx = findall(i -> tlo <= t[i] <= thi, eachindex(t))
    length(idx) >= 12 || return nothing
    tw = Float64.(t[idx]); yw = Float64.(y[idx]); tau = tw .- tw[1]
    extra = tail === nothing ? Vector{Float64}[] : [Float64.(tail[idx])]
    design(G, w) = hcat(exp.(-G .* tau) .* cos.(w .* tau),
                        exp.(-G .* tau) .* sin.(w .* tau), extra...)
    function residual(G, w)
        X = design(G, w)
        c = lstsq_small(X, yw)
        r = copy(yw)
        for j in eachindex(c)
            r .-= c[j] .* view(X, :, j)
        end
        return sum(abs2, r), c
    end
    best = (Inf, 0.0, 0.0)
    Gs = range(0.03/mass, 0.25/mass, length=89)
    ws = range(0.03/mass, 0.25/mass, length=89)
    for pass in 1:3
        for G in Gs, w in ws
            r, _ = residual(G, w)
            r < best[1] && (best = (r, G, w))
        end
        pass == 3 && break
        dG, dw = step(Gs), step(ws)
        Gs = range(max(best[2] - 2dG, 1e-6), best[2] + 2dG, length=41)
        ws = range(max(best[3] - 2dw, 1e-6), best[3] + 2dw, length=41)
    end
    isfinite(best[1]) || return nothing
    rss, c = residual(best[2], best[3])
    return (omegaM = mass*best[3], gammaM = mass*best[2], a = c[1], b = c[2],
            kappa = tail === nothing ? NaN : c[3],
            rel_rms = sqrt(rss/max(sum(abs2, yw), eps(Float64))),
            tlo = tw[1], thi = tw[end], n = length(tw), mass = Float64(mass))
end

end # module Analysis
