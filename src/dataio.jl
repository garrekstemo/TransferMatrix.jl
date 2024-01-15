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

"""
    printstruct(s::Structure, unit=1e9)

Print each layer and its thickness in a somewhat 
visually useful way. Change the default unit multiplier to switch
from nanometers to micrometers. This does not affect any calculations,
only what is printed to the command line when using `printstruct`.
"""
function printstruct(s::Structure, unit=1e9)

    unitstring = "nm"
    if unit == 1e6
        unitstring = "μm"
    end

    print("\n")
    for layer in s.layers

        print("-"^30, "\n")
        print("    $(layer.material), d = $(layer.thickness * unit) $(unitstring)\n")

    end
    print("-"^30, "\n")
end