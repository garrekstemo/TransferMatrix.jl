"""
    load_refractivedata(material::RefractiveMaterial, thickness, wavelength_unit=1e-6)

Retrieves refractive index data from refractiveindex.info database via 
the RefractiveIndex.jl interface.
Wavelengths are in units of micrometers by default.
and two refractive index columns: one for the real part and
the other for the imaginary part.

"""
function load_refractive_data(material::RefractiveMaterial, thickness, wavelength_unit=1.0)
    if haskey(material.data, :data)
        data = readdlm(IOBuffer(material.data[:data]))
        λs = data[:, 1] * wavelength_unit
    else
        println("No data found for $(material.name). Calculate using dispersion formula.")
    end
    if size(data)[2] == 3
        layer = Layer(material.name, thickness, λs, data[:, 2], data[:, 3])
    elseif size(data)[2] == 2
        layer = Layer(material.name, thickness, λs, data[:, 2], zeros(length(data[:, 1])))
        println("No extinction coefficient data for $(material.name). Setting to zero.")
    else
        println("Incorrect data format for $(material.name)")
    end
    return layer
end