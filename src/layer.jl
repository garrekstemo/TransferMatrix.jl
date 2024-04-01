"""
    Layer(material, thickness)

Construct a single layer with keywords:

* `material`: refractive material containing dispersion and extinction data (if available)
* `thickness`: thickness of the layer
"""
struct Layer
    dispersion::Function
    thickness::Real

    function Layer(material, thickness)
        thickness ≥ 0 || throw(DomainError("Layer thickness must be non-negative"))
        new(material, thickness)
    end
end

Layer(material::RefractiveMaterial, thickness::Real) = Layer(refractive_index(material), thickness)
Layer(λs::Vector{Float64}, dispersion::Vector{Real}, extinction::Vector{Real}, thickness::Real) = Layer(refractive_index(λs, dispersion, extinction), thickness)
Layer(dispersion::Function, thickness::Real) = Layer(dispersion, thickness)

"""
    refractive_index()

Return a function that takes a wavelength and gives the real and imaginary parts of the refractive index
"""
function refractive_index(material::RefractiveMaterial)
    return λ -> begin
        try
            return dispersion(material, λ) + im * extinction(material, λ)
        catch e
            if isa(e, ArgumentError)
                return 1.0 + im * 0.0
            else
                rethrow(e)
            end
        end
    end
end
function refractive_index(λs, ns, ks)
    n_real = LinearInterpolation(ns, λs)
    n_imag = LinearInterpolation(ks, λs)
    return λ -> begin
        n_real(λ) + im * n_imag(λ)
    end
end

"""
    find_layerbounds(layers)

Find the unitful z coordinate for all layer-layer interfaces in the structure,
with the first interface starting at z = 0.
(negative z corresponds to positions inside the first layer.)
"""
function find_bounds(layers)

    total_thickness = 0.0
    interface_positions = Float64[]
    
    for layer in layers
        push!(interface_positions, total_thickness + layer.thickness)
        total_thickness += layer.thickness
    end

    return interface_positions, total_thickness
end


"""
    get_refractive_index(material, λ)

Retrieve the refractive index for a material at a given wavelength.
"""
# function get_refractive_index(layer::Layer, λ)
#     n_real = dispersion(layer.material, λ)
#     n_imag = 0.0
    
#     try
#         n_imag = extinction(layer.material, λ)
#     catch e
#         if isa(e, ArgumentError)
#             # Handle the specific case of ArgumentError, which indicates missing extinction data
#             n_imag = 0.0
#         else
#             # Re-throw the error if it's not an ArgumentError
#             rethrow(e)
#         end
#     end
    
#     return n_real + n_imag * im
# end

function get_refractive_index(layer::Layer, λ)
    return layer.dispersion(λ)
end

"""
    dielectric_constant(n_re::Real, n_im::Real)

Return the complex dielectric function from
the real and imaginary parts of the index of refraction.

The complex index of refraction, given by

        n' = n_re + i * n_im
    
(in terms of n_re and n_im), can be used to
obtain the frequency-dependent complex dielectric function

        ε_r(ω) = ε' + iε''

via the relation

        (n_re + i * n_im)^2 = ε' + iε''.
"""
dielectric_constant(n_re::Real, n_im::Real) = (n_re + n_im * im)^2
dielectric_constant(n::Complex) = n^2