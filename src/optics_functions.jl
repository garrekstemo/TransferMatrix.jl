"""
    fresnel(θ, n1, n2)

Calculate the reflectance for s-polarized and p-polarized light
given the incidence angle `θ` (in radians) and indices of refraction
of two media `n1` and `n2` at a plane interface.

Returns `(Rs, Rp)` where `Rs` is s-polarized reflectance and `Rp` is p-polarized reflectance.

The Fresnel equations for reflectance are:

```math
R_s = \\left| \\frac{n_1 \\cos\\theta_i - n_2 \\cos\\theta_t}{n_1 \\cos\\theta_i + n_2 \\cos\\theta_t} \\right|^2
```

```math
R_p = \\left| \\frac{n_2 \\cos\\theta_i - n_1 \\cos\\theta_t}{n_2 \\cos\\theta_i + n_1 \\cos\\theta_t} \\right|^2
```

where ``\\theta_t`` is the transmitted angle given by Snell's law:
``n_1 \\sin\\theta_i = n_2 \\sin\\theta_t``.

Special cases:
- At grazing incidence (θ → π/2), returns `(1.0, 1.0)`
- For total internal reflection (when `n1 > n2` and `θ > θ_critical`), returns `(1.0, 1.0)`

See also: [`fresnel_coefficients`](@ref) for complex amplitude coefficients.

[Fresnel equations](https://en.wikipedia.org/wiki/Fresnel_equations)
"""
function fresnel(θ, n1, n2)
    rs, rp, _, _ = fresnel_coefficients(θ, n1, n2)
    return abs2(rs), abs2(rp)
end

"""
    fresnel_coefficients(θ, n1, n2)

Calculate the Fresnel reflection and transmission amplitude coefficients
for s-polarized and p-polarized light at a plane interface.

Returns `(rs, rp, ts, tp)` where:
- `rs`, `rp`: reflection amplitude coefficients (complex)
- `ts`, `tp`: transmission amplitude coefficients (complex)

These are the amplitude (not intensity) coefficients, so reflectance `R = |r|²`
and transmittance requires the full expression `T = (n2 cosθt)/(n1 cosθi) |t|²`.

For total internal reflection, the reflection coefficients have unit magnitude
with a phase shift (evanescent wave), and transmission coefficients are zero.

See also: [`fresnel`](@ref) for intensity reflectances.
"""
function fresnel_coefficients(θ, n1, n2)
    cosθ = cos(θ)
    sinθ = sin(θ)

    # Handle grazing incidence
    if abs(cosθ) < eps(Float64)
        # At grazing: total reflection, no transmission
        return (complex(1.0), complex(1.0), complex(0.0), complex(0.0))
    end

    sin_ratio = n1 * sinθ / n2
    sin_ratio_sq = sin_ratio^2

    if real(sin_ratio_sq) > 1.0 && imag(n2) ≈ 0
        # Total internal reflection: reflection with phase shift
        # cosθt becomes purely imaginary
        cosθt = im * √(sin_ratio_sq - 1)

        n1_cosθ = n1 * cosθ
        n2_cosθt = n2 * cosθt
        n2_cosθ = n2 * cosθ
        n1_cosθt = n1 * cosθt

        rs = (n1_cosθ - n2_cosθt) / (n1_cosθ + n2_cosθt)
        rp = (n2_cosθ - n1_cosθt) / (n2_cosθ + n1_cosθt)

        # Evanescent transmission (effectively zero power transfer)
        return (rs, rp, complex(0.0), complex(0.0))
    end

    cosθt = √(complex(1 - sin_ratio_sq))

    n1_cosθ = n1 * cosθ
    n2_cosθt = n2 * cosθt
    n2_cosθ = n2 * cosθ
    n1_cosθt = n1 * cosθt

    # Reflection amplitudes
    rs = (n1_cosθ - n2_cosθt) / (n1_cosθ + n2_cosθt)
    rp = (n2_cosθ - n1_cosθt) / (n2_cosθ + n1_cosθt)

    # Transmission amplitudes
    ts = 2 * n1_cosθ / (n1_cosθ + n2_cosθt)
    tp = 2 * n1_cosθ / (n2_cosθ + n1_cosθt)

    return (rs, rp, ts, tp)
end

"""
    airy(n0, nf, ns, d, λ; θ=0.0)

Calculate the reflectance and transmittance of a single thin film using the
exact Airy formula (multiple-beam interference).

# Arguments
- `n0`: refractive index of incident medium (can be complex)
- `nf`: refractive index of the film (can be complex)
- `ns`: refractive index of the substrate (can be complex)
- `d`: film thickness (same units as λ)
- `λ`: wavelength (same units as d)
- `θ=0.0`: angle of incidence in radians

# Returns
`(Rs, Rp, Ts, Tp)` - reflectance and transmittance for s and p polarizations.

# Physics
The Airy formula accounts for all multiple reflections within the film:

```math
r = \\frac{r_{01} + r_{12} e^{2i\\delta}}{1 + r_{01} r_{12} e^{2i\\delta}}
```

where `δ = 2π nf d cos(θf) / λ` is the phase thickness and `r₀₁`, `r₁₂` are
the Fresnel reflection coefficients at the two interfaces.

This provides an exact analytical solution for validating numerical TMM results.

# Example
```julia
# Quarter-wave anti-reflection coating
n_air, n_film, n_glass = 1.0, 1.38, 1.52
λ = 0.55  # μm
d = λ / (4 * n_film)  # quarter-wave thickness
Rs, Rp, Ts, Tp = airy(n_air, n_film, n_glass, d, λ)
```

See also: [`fresnel`](@ref), [`fresnel_coefficients`](@ref)
"""
function airy(n0, nf, ns, d, λ; θ=0.0)
    sinθ = sin(θ)
    cosθ = cos(θ)

    # Snell's law for angles in film and substrate
    sinθf = n0 * sinθ / nf
    sinθs = n0 * sinθ / ns

    # Handle potential TIR or complex angles
    cosθf = √(complex(1 - sinθf^2))
    cosθs = √(complex(1 - sinθs^2))

    # Phase thickness in the film
    δ = 2π * nf * d * cosθf / λ

    # Fresnel amplitudes at interface 0→f (incident → film)
    rs_01, rp_01, ts_01, tp_01 = fresnel_coefficients(θ, n0, nf)

    # Fresnel amplitudes at interface f→s (film → substrate)
    # Need angle in film for this interface
    θf = asin(sinθf)
    rs_12, rp_12, ts_12, tp_12 = fresnel_coefficients(real(θf), nf, ns)

    # Airy formula for total reflection amplitude
    # r = (r01 + r12 * exp(2iδ)) / (1 + r01 * r12 * exp(2iδ))
    phase = exp(2im * δ)

    rs_total = (rs_01 + rs_12 * phase) / (1 + rs_01 * rs_12 * phase)
    rp_total = (rp_01 + rp_12 * phase) / (1 + rp_01 * rp_12 * phase)

    # Airy formula for total transmission amplitude
    # t = (t01 * t12 * exp(iδ)) / (1 + r01 * r12 * exp(2iδ))
    phase_single = exp(im * δ)

    ts_total = (ts_01 * ts_12 * phase_single) / (1 + rs_01 * rs_12 * phase)
    tp_total = (tp_01 * tp_12 * phase_single) / (1 + rp_01 * rp_12 * phase)

    # Reflectance is simply |r|²
    Rs = abs2(rs_total)
    Rp = abs2(rp_total)

    # Transmittance requires correction for different media
    # T = (ns * cosθs) / (n0 * cosθ) * |t|²
    # For real indices, this ensures energy conservation: R + T = 1
    prefactor_s = real(ns * cosθs) / real(n0 * cosθ)
    prefactor_p = real(ns * cosθs) / real(n0 * cosθ)

    Ts = prefactor_s * abs2(ts_total)
    Tp = prefactor_p * abs2(tp_total)

    return (Rs, Rp, Ts, Tp)
end

"""
    stopband(n1, n2)

Calculate the frequency bandwidth Δf of the photonic stopband
for a distributed bragg reflector (DBR) with two alternating
materials of refractive indices `n1` and `n2`.

```math
    \\frac{\\Delta f_0}{f_0} = \\frac{4}{\\pi} \\arcsin \\left( \\frac{n_2 - n_1}{n_2 + n_1} \\right)
```

[Distributed Bragg reflector](https://en.wikipedia.org/wiki/Distributed_Bragg_reflector)
"""
function stopband(n1, n2)
    4 * asin(abs(n2 - n1) / (n2 + n1)) / π
end

"""
    dbr_reflectivity(no, ns, n1, n2, N)

Approximate the reflectivity of a DBR structure with originating medium with refractive index `no`,
substrate with index `ns`, and alternating materials with indices `n1` and `n2` and number of repetitions
`N`. The repeated pair of materials are assumed to have quarter-wave thickness ``nd = \\lambda / 4``,
where ``n`` is the refractive index, ``d`` is the layer thickness, and ``\\lambda`` is the wavelength of the light.

[Distributed Bragg reflector](https://en.wikipedia.org/wiki/Distributed_Bragg_reflector)
"""
function dbr_reflectivity(no, ns, n1, n2, N)
    r = (no * n2^(2 * N) - ns * n1^(2 * N)) / (no * n2^(2 * N) + ns * n1^(2 * N))
    return r^2
end