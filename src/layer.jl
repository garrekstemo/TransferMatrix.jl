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
Layer(λs::AbstractVector, dispersion::AbstractVector, extinction::AbstractVector, thickness::Real) = Layer(refractive_index(λs, dispersion, extinction), thickness)

"""
    refractive_index()

Return a function that takes a wavelength and gives the real and imaginary parts of the refractive index
"""
function refractive_index(material::RefractiveMaterial)
    return λ -> begin
        n_imag = 0.0im  # Define n_imag before the try block
        try
            n_imag = im * extinction(material, λ)
        catch e
            if isa(e, ArgumentError)
                n_imag = 0.0im
            else
                rethrow(e)
            end
        end
        return dispersion(material, λ) + n_imag
    end
end
function refractive_index(λs::AbstractVector, ns::AbstractVector, ks::AbstractVector)
    n_real = LinearInterpolation(ns, λs)
    n_imag = LinearInterpolation(ks, λs)
    return λ -> begin
        n_real(λ) + im * n_imag(λ)
    end
end

"""
    find_bounds(layers)

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