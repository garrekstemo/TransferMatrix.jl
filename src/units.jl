"""
    _to_um(x)

Convert a length quantity to micrometers (the package's internal length unit).

No-op for plain `Real` inputs (already assumed to be μm). The `UnitfulExt`
extension adds a method for `Unitful.Length` that strips to μm. Used to normalize
layer thickness and the field-sampling step `dz`.
"""
_to_um(x) = x

"""
    _to_wavelength_um(x)

Convert a spectral input to wavelength in micrometers.

No-op for plain `Real` inputs (already assumed to be μm). The `UnitfulExt`
extension adds a method that accepts a `Unitful` length, spectroscopic wavenumber
(ν̃ = 1/λ), ordinary frequency (λ = c/f), or photon energy (λ = hc/E). Used to
normalize the wavelength `λ` at every public entry point.
"""
_to_wavelength_um(x) = x

"""
    _to_radians(x)

Convert an angle to radians (the package's internal angle unit).

No-op for plain `Real` inputs (already assumed to be radians). The `UnitfulExt`
extension adds a method for unit-bearing angles (e.g. `45u"°"`) that strips to
radians. Used to normalize the incidence angle `θ` at every public entry point.
"""
_to_radians(x) = x

"""
    _to_eV(x)

Convert a dispersion-model frequency/energy parameter (`ω_p`, `ω_0`, `γ`) to
electron-volts (the internal unit for the built-in dispersion models).

No-op for plain `Real` inputs (already assumed to be eV). The `UnitfulExt`
extension adds methods for unit-bearing energy (`eV`, `J`), spectroscopic
wavenumber (`cm^-1`), and ordinary frequency (`Hz`), each mapped *linearly* to
energy so resonance positions and damping widths convert identically.
"""
_to_eV(x) = x
