# Analytic dispersion models (Drude metals, Lorentz oscillators).
#
# Each public function returns a closure `λ -> n` giving the complex refractive
# index at vacuum wavelength `λ` (μm). The permittivity is evaluated in
# frequency space and the μm wavelength is converted to photon energy (eV)
# internally, so all frequency-like parameters share the eV convention.
#
# Time convention: exp(-iωt). Signs are chosen so ε″ > 0 (absorption).

"""
    HC_EV_UM

Product `hc/e` in eV·μm (CODATA 2018). Converts a vacuum wavelength in μm to
photon energy in eV via `E = HC_EV_UM / λ`.
"""
const HC_EV_UM = 6.62607015e-34 * c_0 / 1.602176634e-19 * 1e6

# Vacuum wavelength (μm) → photon energy (eV).
_photon_energy_eV(λ) = HC_EV_UM / λ

# Drude permittivity contribution (relative to ε∞): -ω_p² / (ω² + iγω).
_drude_eps(ωp, γ, ω) = -ωp^2 / (ω^2 + im * γ * ω)

"""
    drude(ω_p, γ; ε_inf=1.0)

Return a closure `λ -> n` for the Drude free-carrier (metal) model:

```math
ε(ω) = ε_\\infty - \\frac{ω_p^2}{ω^2 + iγω}, \\qquad n = \\sqrt{ε}.
```

# Arguments
- `ω_p`: plasma energy (eV by default; a `Unitful` energy/wavenumber/frequency is
  normalized to eV when `Unitful` is loaded).
- `γ`: damping/collision rate (eV).
- `ε_inf`: high-frequency dielectric constant (dimensionless, default `1.0`).

The returned closure takes a vacuum wavelength `λ` in **μm** and returns the
complex refractive index `n` (with `n″ > 0` for absorption under the package's
`exp(-iωt)` convention).

# Example
```julia
# Gold-like Drude metal, 50 nm film
au = drude(9.0, 0.07)        # ω_p = 9 eV, γ = 0.07 eV
layer = Layer(au, 0.05)
```

See also: [`lorentz`](@ref), [`drude_lorentz`](@ref).
"""
function drude(ω_p, γ; ε_inf=1.0)
    ωp = _to_eV(ω_p)
    g = _to_eV(γ)
    return λ -> begin
        ω = _photon_energy_eV(λ)
        return sqrt(ε_inf + _drude_eps(ωp, g, ω))
    end
end

# Summed Lorentz-oscillator permittivity contribution (relative to ε∞). Each
# term is (ω0, Δε, γ) with ω0, γ in eV and Δε dimensionless.
function _lorentz_eps(terms, ω)
    χ = 0.0im
    for (ω0, Δε, γ) in terms
        χ += Δε * ω0^2 / (ω0^2 - ω^2 - im * γ * ω)
    end
    return χ
end

"""
    lorentz(ω_0, Δε, γ; ε_inf=1.0)
    lorentz(oscillators; ε_inf=1.0)

Return a closure `λ -> n` for the Lorentz oscillator model:

```math
ε(ω) = ε_\\infty + \\sum_j \\frac{Δε_j\\, ω_{0j}^2}{ω_{0j}^2 - ω^2 - iγ_j ω},
\\qquad n = \\sqrt{ε}.
```

# Arguments
- `ω_0`: resonance energy (eV by default).
- `Δε`: dimensionless oscillator strength (the static contribution of the mode;
  `ε(0) = ε_inf + Σ Δε_j`).
- `γ`: damping (eV).
- `ε_inf`: high-frequency dielectric constant (default `1.0`).
- `oscillators`: an iterable of `(ω_0, Δε, γ)` 3-tuples for a multi-oscillator
  model (e.g. several phonon modes).

The returned closure takes a vacuum wavelength `λ` in **μm**.

# Examples
```julia
# Single oscillator
n = lorentz(2.0, 3.0, 0.05)                       # ω_0 = 2 eV

# Two oscillators on top of ε∞ = 2.5
n = lorentz([(2.0, 1.0, 0.05), (3.5, 0.4, 0.1)]; ε_inf=2.5)
```

See also: [`drude`](@ref), [`drude_lorentz`](@ref).
"""
function lorentz(oscillators; ε_inf=1.0)
    terms = [(_to_eV(ω0), Δε, _to_eV(γ)) for (ω0, Δε, γ) in oscillators]
    return λ -> begin
        ω = _photon_energy_eV(λ)
        return sqrt(ε_inf + _lorentz_eps(terms, ω))
    end
end

lorentz(ω_0, Δε, γ; ε_inf=1.0) = lorentz(((ω_0, Δε, γ),); ε_inf=ε_inf)

"""
    drude_lorentz(ω_p, γ_D, oscillators; ε_inf=1.0)

Return a closure `λ -> n` for a combined Drude + Lorentz model (the standard form
for real metals such as gold):

```math
ε(ω) = ε_\\infty - \\frac{ω_p^2}{ω^2 + iγ_D ω}
       + \\sum_j \\frac{Δε_j\\, ω_{0j}^2}{ω_{0j}^2 - ω^2 - iγ_j ω}.
```

# Arguments
- `ω_p`, `γ_D`: Drude plasma energy and damping (eV).
- `oscillators`: iterable of Lorentz `(ω_0, Δε, γ)` 3-tuples (may be empty for a
  pure Drude metal).
- `ε_inf`: high-frequency dielectric constant (default `1.0`).

The returned closure takes a vacuum wavelength `λ` in **μm**.

See also: [`drude`](@ref), [`lorentz`](@ref).
"""
function drude_lorentz(ω_p, γ_D, oscillators; ε_inf=1.0)
    ωp = _to_eV(ω_p)
    gD = _to_eV(γ_D)
    terms = [(_to_eV(ω0), Δε, _to_eV(γ)) for (ω0, Δε, γ) in oscillators]
    return λ -> begin
        ω = _photon_energy_eV(λ)
        return sqrt(ε_inf + _drude_eps(ωp, gD, ω) + _lorentz_eps(terms, ω))
    end
end
