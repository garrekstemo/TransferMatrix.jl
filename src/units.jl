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
