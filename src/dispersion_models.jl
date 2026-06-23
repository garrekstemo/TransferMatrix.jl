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
