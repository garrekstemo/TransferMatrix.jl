"""
    read_refractive(f::String, material::String, thickness; div=1.0, freq=false)

Reads a csv file from the database website refractiveindex.info containing 
a wavelength column in units of micrometers,
and two refractive index columns: one for the real part and
the other for the imaginary part.

The header names are:

"Wavelength, μm", "n", "k"

Note that if you put a file in the `refractive_index_data` folder, then
this function can automatically find it if you simply put the file name
(with .csv extension) in the `refractive_filename` field in the yaml file.
If you put the file somewhere else, you must provide the full path to that file
(or relative path).

For some reason, the wavelength column cannot be normalized (so as to be a valid Julia identifier)
using `normalizednames = true` in `CSV.File`, so the header is skipped.
"""
function read_refractive(f::String, material::String, thickness; div=1.0, freq=false)
    if !(isfile(f))
        f = abspath("./refractive_index_data/$(f)")
    else
        f = abspath(f)
    end
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