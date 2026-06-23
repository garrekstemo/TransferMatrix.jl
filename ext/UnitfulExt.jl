module UnitfulExt

using TransferMatrix
using Unitful
import TransferMatrix: _to_um, _to_wavelength_um, _to_radians, _to_eV

# Length (thickness, dz) → μm as Float64.
_to_um(x::Unitful.Length) = Float64(ustrip(u"μm", x))

# Angle of incidence (θ) → radians as Float64. Angles are dimensionless in
# Unitful, so this catches both u"°" and u"rad"; a plain Real falls through to
# the core no-op and is assumed to already be in radians.
_to_radians(x::Unitful.DimensionlessQuantity) = Float64(ustrip(u"rad", x))

# Spectral input (λ) → wavelength in μm as Float64. Dispatch by physical
# dimension. Only the length dimension is handled here; spectral dimensions are
# added in a later step.
function _to_wavelength_um(x::Unitful.Quantity)
    d = dimension(x)
    d == dimension(u"m")    && return Float64(ustrip(u"μm", x))
    d == dimension(u"m^-1") && return Float64(ustrip(u"μm", 1 / x))
    d == dimension(u"Hz")   && return Float64(ustrip(u"μm", Unitful.c0 / x))
    d == dimension(u"eV")   && return Float64(ustrip(u"μm", Unitful.h * Unitful.c0 / x))
    throw(ArgumentError("cannot interpret $x as a wavelength or spectral quantity"))
end

# Dispersion-model parameter (ω_p, ω_0, γ) → energy in eV as Float64. Dispatch by
# physical dimension with the *linear* frequency→energy map (the mirror of the
# inverse _to_wavelength_um), so resonance positions and damping widths convert
# identically. Reuses Unitful.h / Unitful.c0.
function _to_eV(x::Unitful.Quantity)
    d = dimension(x)
    d == dimension(u"eV")   && return Float64(ustrip(u"eV", x))
    d == dimension(u"m^-1") && return Float64(ustrip(u"eV", Unitful.h * Unitful.c0 * x))
    d == dimension(u"Hz")   && return Float64(ustrip(u"eV", Unitful.h * x))
    throw(ArgumentError("cannot interpret $x as a dispersion-model frequency/energy parameter"))
end

# Construct a Layer with a unit-bearing thickness; strip to μm so the stored
# field stays a plain Float64. Covers three forms:
#   - isotropic 2-arg: Layer(material, t)
#   - anisotropic 3-dispersion: Layer(nx, ny, nz, t)
#   - tabulated-data: Layer(λs, ns, ks, t)
# The tabulated method (all-AbstractVector) is strictly more specific than the
# anisotropic (Any, Any, Any, Length) method, so it is selected first with no
# ambiguity when all three positional args are AbstractVectors.
TransferMatrix.Layer(material, t::Unitful.Length) = TransferMatrix.Layer(material, _to_um(t))
TransferMatrix.Layer(nx, ny, nz, t::Unitful.Length; euler=(0.0, 0.0, 0.0)) =
    TransferMatrix.Layer(nx, ny, nz, _to_um(t); euler=euler)
TransferMatrix.Layer(λs::AbstractVector, ns::AbstractVector, ks::AbstractVector, t::Unitful.Length) =
    TransferMatrix.Layer(λs, ns, ks, _to_um(t))

end # module
