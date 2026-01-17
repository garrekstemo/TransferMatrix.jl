"""
    Layer(material, thickness)

Construct a single layer for transfer matrix calculations.

# Arguments
- `material`: Refractive material (from RefractiveIndex.jl) or a dispersion function `λ -> n(λ)`
- `thickness`: Layer thickness in the same units as wavelength (typically μm)

# Units Convention
All length quantities (wavelength `λ`, thickness, `dz`) must use consistent units.
Following RefractiveIndex.jl conventions, **micrometers (μm)** are recommended.

`Layer` is parametric as `Layer{F,T}` where `F` is the dispersion function type and `T` is the thickness type.

# Example
```julia
# Using RefractiveIndex.jl material (wavelength in μm)
n_sio2 = RefractiveMaterial("main", "SiO2", "Malitson")
layer = Layer(n_sio2, 0.1)  # 100 nm = 0.1 μm

# Using custom dispersion function
layer = Layer(λ -> 1.5 + 0.01im, 0.25)  # constant n = 1.5 + 0.01i
```
"""
struct Layer{F,T<:Real}
    dispersion::F
    thickness::T

    function Layer(material::F, thickness::T) where {F,T<:Real}
        thickness ≥ 0 || throw(DomainError("Layer thickness must be non-negative"))
        new{F,T}(material, thickness)
    end
end

Layer(material::RefractiveMaterial, thickness::Real) = Layer(refractive_index(material), thickness)
Layer(λs::AbstractVector, dispersion::AbstractVector, extinction::AbstractVector, thickness::Real) = Layer(refractive_index(λs, dispersion, extinction), thickness)

"""
    refractive_index(material::RefractiveMaterial)

Return a function that takes a wavelength and gives the complex refractive index.

The extinction coefficient availability is checked once at construction time to avoid
try-catch overhead in the inner loop of spectral calculations.
"""
function refractive_index(material::RefractiveMaterial)
    # Check at construction time whether extinction data is available
    has_extinction = try
        extinction(material, 1.0)  # Probe with test wavelength
        true
    catch e
        isa(e, ArgumentError) ? false : rethrow(e)
    end

    if has_extinction
        return λ -> dispersion(material, λ) + im * extinction(material, λ)
    else
        return λ -> dispersion(material, λ) + 0.0im
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
    n_layers = length(layers)
    interface_positions = Vector{Float64}(undef, n_layers)

    total_thickness = 0.0
    for (i, layer) in enumerate(layers)
        total_thickness += layer.thickness
        interface_positions[i] = total_thickness
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
