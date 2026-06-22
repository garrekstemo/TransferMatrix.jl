module RefractiveIndexExt

using TransferMatrix
using RefractiveIndex
import TransferMatrix: refractive_index

# Build a complex dispersion function n(λ) = n(λ) + i·k(λ) from a RefractiveIndex.jl
# material. Extinction-data availability is probed once here, not in the spectral
# inner loop: an ArgumentError means the material has no k data; any other error
# means k exists but the probe wavelength was out of range.
function refractive_index(material::RefractiveMaterial)
    has_extinction = try
        RefractiveIndex.extinction(material, 1.0)
        true
    catch e
        e isa ArgumentError ? false : true
    end
    if has_extinction
        return λ -> RefractiveIndex.dispersion(material, λ) + im * RefractiveIndex.extinction(material, λ)
    else
        return λ -> RefractiveIndex.dispersion(material, λ) + 0.0im
    end
end

# Convenience constructors that funnel a material through `refractive_index`. Each
# extends a TransferMatrix-owned function (no type piracy) and is strictly more
# specific than its core fallback, so it is selected without ambiguity.
TransferMatrix.Layer(material::RefractiveMaterial, thickness::Real) =
    TransferMatrix.Layer(refractive_index(material), thickness)

TransferMatrix._index_fn(m::RefractiveMaterial) = refractive_index(m)

TransferMatrix.Sheet(material::RefractiveMaterial, d::Real) =
    TransferMatrix._diagonal_sheet(TransferMatrix._sigma_from_index(refractive_index(material), d))

end # module
