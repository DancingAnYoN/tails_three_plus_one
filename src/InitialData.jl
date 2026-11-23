# Constraint-satisfying initial slice: Schwarzschild black hole of mass M0 plus
# an ingoing Gaussian scalar pulse, at t = 0 with f(0) = 0 (so J = r).
module InitialData
using ..Types: Config, IA, IB, IPHI, IPI, IFIELD
using ..Grid: radial_grid
using ..Geometry: schwarzschild_fields, schwarzschild_Bprime
using ..ScalarField: gaussian
export initial_data
"""
Initial constraint ODEs for delta_m=m-M0 and delta_B=B-B_Schw(M0),
with the Marsa--Choptuik Schwarzschild-inspired A(m,r) ansatz. Evolving
only differences exactly preserves the vacuum initial-data background.
Ingoing scalar data use ∂t(r φ)=∂r(r φ), f(0)=fdot(0)=0.

Returns (∂_r δm, ∂_r δB): the Hamiltonian constraint in mass form gives ∂_r m,
and the momentum constraint gives ∂_r B once A is fixed by the ansatz
A = -2m(r+m)/[r(r+2m)]^{3/2} (Schwarzschild A with M replaced by the local m).
"""
function initial_ode(r::Float64, dm::Float64, dB::Float64, p::Config)
    m = p.M0 + dm
    B = schwarzschild_fields(r,p.M0).B + dB
    inva2 = 1 - 2m/r + (r*B)^2
    (isfinite(inva2) && inva2 > 0) || error("Non-positive initial inverse a² at r=$r")
    a = inv(sqrt(inva2))
    A = -2m*(r+m)/(r*(r+2m))^1.5
    phi, Phi = gaussian(r,p)
    q = a*r*B
    (1+q > 0) || error("Initial gauge denominator failed")
    # Ingoing wave r φ = g(t + r), since v = t + r is the ingoing null coordinate:
    # ∂_t φ = Φ + φ/r, hence Π = (∂_t φ - βΦ)/(1-β) = Φ + φ/(r(1-β)).
    # Pi = Phi + phi/(r(1-beta)); 1/(1-beta) = 1+q at fdot=0.
    Pi = Phi + (1+q)*phi/r
    mp = 4π*r*r*((Phi*Phi+Pi*Pi)/(2a*a) + (r*B/a)*Phi*Pi)
    Bp = (A-B)/r + 4π*Phi*Pi/a - schwarzschild_Bprime(r,p.M0)
    return mp,Bp
end

"""
Build the initial state `(r, u)` on the grid of `p`.

Integrates (δm, δB) outward from the horizon, where δm = δB = 0 so that
m(r_h) = M0, B(r_h) = 1/(2√2 M0) and aJB = 1. The scalar field vanishes inside
the cutoff, so for a vacuum run (ε = 0) the slice is exactly Schwarzschild.
"""
function initial_data(p::Config)
    r = radial_grid(p)
    n = length(r)
    u = Matrix{Float64}(undef, n, 5)
    dm = Vector{Float64}(undef, n)
    dB = Vector{Float64}(undef, n)
    dm[1] = 0.0
    dB[1] = 0.0 # a_h J_h B_h=1 for the Schwarzschild reference
    h = p.dr
    # Classical RK4 in r for the two coupled constraint ODEs.
    @inbounds for i in 1:n-1
        x = r[i]
        m0, b0 = dm[i], dB[i]
        k1m,k1b = initial_ode(x,m0,b0,p)
        k2m,k2b = initial_ode(x+h/2,m0+h*k1m/2,b0+h*k1b/2,p)
        k3m,k3b = initial_ode(x+h/2,m0+h*k2m/2,b0+h*k2b/2,p)
        k4m,k4b = initial_ode(x+h,m0+h*k3m,b0+h*k3b,p)
        # `2k2m` means 2*k2m (identifier, not a numeric literal).
        dm[i+1] = m0 + h*(k1m+2k2m+2k3m+k4m)/6
        dB[i+1] = b0 + h*(k1b+2k2b+2k3b+k4b)/6
    end
    # Assemble the state columns from m, B and the pulse.
    @inbounds for i in 1:n
        J = r[i]
        mi = p.M0 + dm[i]
        Bi = schwarzschild_fields(J,p.M0).B + dB[i]
        a = inv(sqrt(1 - 2mi/J + (J*Bi)^2))
        phi,Phi = gaussian(J,p)
        Pi = Phi + (1+a*J*Bi)*phi/J
        u[i,IA] = a
        u[i,IB] = Bi
        u[i,IPHI] = Phi
        u[i,IPI] = Pi
        u[i,IFIELD] = phi
    end
    return r,u
end

end # module InitialData
