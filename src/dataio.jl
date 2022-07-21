const ε_0 = 8.8541878128e-12
const μ_0 = 1.25663706212e-6
const c_0 = 299792458

"""
Reads a csv file from refractiveindex.info containing 
a wavlength column in units of micrometers,
a real refractive index column, and an imaginary refractive index column.
The header names are:

"Wavelength, μm", "n", "k"

For some reason, the wavelength column cannot be normalized (so as to be a valid Julia identifier)
using `normalizednames = true` in `CSV.File`, so the header is skipped.
"""
function read_refractive(f::String, material::String, thickness::Float64; div::Float64=1.0, freq::Bool=false)

    ndata = CSV.File(f, skipto = 2, header = false, types = Float64)
    if freq == true
        λ = c_0 ./ ndata.Column1
    else
        λ = ndata.Column1 ./ div
    end

    if length(ndata.names) == 3
        layer = Layer(material, thickness, λ, ndata.Column2, ndata.Column3)
    elseif length(ndata.names) == 2
        layer = Layer(material, thickness, λ, ndata.Column2, fill(0.0, length(λ)))
    end
    return layer
end
