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
# field stays a plain Float64. Covers the isotropic 2-arg and anisotropic 4-arg
# forms. The tabulated-data + unitful-thickness combo
# (Layer(λs, ns, ks, d::Unitful.Length)) is intentionally not wired — use a
# numeric μm thickness there.
TransferMatrix.Layer(material, t::Unitful.Length) = TransferMatrix.Layer(material, _to_um(t))
TransferMatrix.Layer(nx, ny, nz, t::Unitful.Length; euler=(0.0, 0.0, 0.0)) =
    TransferMatrix.Layer(nx, ny, nz, _to_um(t); euler=euler)

end # module
