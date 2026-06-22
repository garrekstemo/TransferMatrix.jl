module UnitfulExt

using TransferMatrix
using Unitful
import TransferMatrix: _to_um, _to_wavelength_um

# Length (thickness, dz) → μm as Float64.
_to_um(x::Unitful.Length) = Float64(ustrip(u"μm", x))

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
