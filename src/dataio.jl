"""
    load_refractivedata(material::RefractiveMaterial, thickness, wavelength_unit=1e-6)

Retrieves refractive index data from refractiveindex.info database via 
the RefractiveIndex.jl interface.
Wavelengths are in units of micrometers by default.
and two refractive index columns: one for the real part and
the other for the imaginary part.

"""
function load_refractive_data(material::RefractiveMaterial, thickness, wavelength_unit=1)
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
    read_yaml(file::String)

Read a yaml file with Structure parameters to
perform transfer matrix calculations.
"""
function load_from_yaml(yamlfile::String, unitconvert = 1.0)

    yml = YAML.load_file(yamlfile)

    λ_min = yml["min_wavelength"] * unitconvert
    λ_max = yml["max_wavelength"] * unitconvert
    n_λ = yml["n_wavelengths"]
    
    θ_min = yml["theta_i"]
    θ_max = yml["theta_f"]
    n_θ = yml["n_theta"]
    θs = range(θ_min, θ_max, length = n_θ)

    s = Structure()
    layers = Layer[]

    # A periodic structure is considered one layer.
    for i in 1:length(yml["layers"])

        layeritem = yml["layers"]["layer$(i)"]
    
        if "periods" in keys(layeritem)
            num_periods = layeritem["periods"]
    
            # Now we deal with the layers within the periodic structure.
            for n in 1:num_periods
    
                for j in 1:length(layeritem["layers"])
    
                    perioditem = layeritem["layers"]["layer$(j)"]
                    material = perioditem["material"]
                    thickness = perioditem["thickness"] * unitconvert
    
                    if "refractive_filename" in keys(perioditem)
                        layer = read_refractive(perioditem["refractive_filename"], material, thickness, div=1e6)
                    
                    elseif "refractive_index" in keys(perioditem)
                        λs, ns, κs = parse_refractive(perioditem)
                        layer = Layer(material, thickness, λs, ns, κs)
    
                    else
                        println("Incorrect yaml config format for periodic layer. See example files.")
    
                    end
                    push!(layers, layer)
                end
            end
        else
            material = layeritem["material"]
            thickness = layeritem["thickness"] * unitconvert

            if "refractive_filename" in keys(layeritem)
                layer = read_refractive(layeritem["refractive_filename"], material, thickness, div=1e6)
    
            elseif "refractive_index" in keys(layeritem)
                λs, ns, κs = parse_refractive(layeritem)
                layer = Layer(material, thickness, λs, ns, κs)
            else
                println("Incorrect yaml config format for layer $(material). See example files.")
            end
    
            push!(layers, layer)
    
        end
        s.layers = layers
    end

    s.layers = layers
    s.θ = θs
    new_λs = range(λ_min, λ_max, length = n_λ)

    return TransferMatrix.initialize(s, new_λs)
end

"""
    yaml_example()

Generate an example input parameter yaml file.
"""
function yaml_example()
end

function parse_refractive(d::Dict)

    λs = d["wavelength"] * 10^-6
    ns = d["refractive_index"]
    κs = d["extinction_coeff"]

    if typeof(λs) != Vector
        λs = Float64[λs]
    end
    if typeof(ns) != Vector
        ns = Float64[ns]
    end
    if typeof(κs) != Vector
        κs = Float64[κs]
    end

    return λs, ns, κs
end