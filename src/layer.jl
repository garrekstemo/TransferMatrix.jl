"""
    Layer(material, thickness)
    Layer(nx, ny, nz, thickness)

Construct a single layer for transfer matrix calculations.

# Isotropic Layer (single refractive index)
```julia
Layer(material, thickness)
```
- `material`: Refractive material (from RefractiveIndex.jl) or a dispersion function `λ -> n(λ)`
- `thickness`: Layer thickness in the same units as wavelength (typically μm)

# Anisotropic Layer (biaxial: three refractive indices)
```julia
Layer(nx, ny, nz, thickness)
```
- `nx`, `ny`, `nz`: Dispersion functions `λ -> n(λ)` for each principal axis
- `thickness`: Layer thickness in μm

For uniaxial materials, set two axes equal (e.g., `nx = ny` for optic axis along z).

# Units Convention
All length quantities (wavelength `λ`, thickness, `dz`) must use consistent units.
Following RefractiveIndex.jl conventions, **micrometers (μm)** are recommended.

`Layer` is parametric as `Layer{F,T}` where `F` is the dispersion function type and `T` is the thickness type.
For anisotropic layers, `F` is a `Tuple` of three dispersion functions.

# Examples
```julia
# Isotropic: Using RefractiveIndex.jl material
n_sio2 = RefractiveMaterial("main", "SiO2", "Malitson")
layer = Layer(n_sio2, 0.1)  # 100 nm = 0.1 μm

# Isotropic: Using custom dispersion function
layer = Layer(λ -> 1.5 + 0.01im, 0.25)  # constant n = 1.5 + 0.01i

# Anisotropic: Uniaxial crystal (calcite-like, optic axis along z)
no = λ -> 1.658  # ordinary index
ne = λ -> 1.486  # extraordinary index
layer = Layer(no, no, ne, 0.5)

# Anisotropic: Biaxial crystal
layer = Layer(λ -> 1.5, λ -> 1.6, λ -> 1.7, 0.3)
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

# Isotropic constructors
Layer(material::RefractiveMaterial, thickness::Real) = Layer(refractive_index(material), thickness)
Layer(λs::AbstractVector, dispersion::AbstractVector, extinction::AbstractVector, thickness::Real) = Layer(refractive_index(λs, dispersion, extinction), thickness)

# Anisotropic constructor: three dispersion functions for nx, ny, nz
"""
    Layer(nx, ny, nz, thickness; euler=(0,0,0))

Construct an anisotropic layer with different refractive indices along each principal axis.

# Arguments
- `nx`, `ny`, `nz`: Dispersion functions `λ -> n(λ)` for each principal axis
- `thickness`: Layer thickness in μm
- `euler`: Optional tuple `(φ, θ, ψ)` of ZYZ Euler angles in radians (default: no rotation)

The dielectric tensor is constructed as a diagonal matrix with ε_i = n_i² in the crystal
frame, then rotated to the lab frame using the Euler angles.

# Euler Angle Convention (ZYZ)
- `φ` (phi): First rotation about z-axis
- `θ` (theta): Rotation about new y-axis (tilt angle from z)
- `ψ` (psi): Second rotation about new z-axis

# Examples
```julia
# Uniaxial crystal with optic axis along z (no rotation needed)
layer = Layer(no, no, ne, 0.5)

# Same crystal with optic axis tilted 45° from z in the xz-plane
layer = Layer(no, no, ne, 0.5; euler=(0, π/4, 0))

# Optic axis in the xy-plane at 30° from x
layer = Layer(no, no, ne, 0.5; euler=(π/6, π/2, 0))
```
"""
function Layer(nx::F1, ny::F2, nz::F3, thickness::T; euler::NTuple{3,Real}=(0.0, 0.0, 0.0)) where {F1,F2,F3,T<:Real}
    thickness ≥ 0 || throw(DomainError("Layer thickness must be non-negative"))
    φ, θ, ψ = Float64.(euler)
    Layer((nx, ny, nz, φ, θ, ψ), thickness)
end

"""
    isanisotropic(layer::Layer)

Return `true` if the layer has anisotropic optical properties (different refractive
indices along principal axes), `false` for isotropic layers.
"""
isanisotropic(layer::Layer) = layer.dispersion isa Tuple

"""
    get_refractive_indices(layer::Layer, λ)

Return the refractive indices for a layer at wavelength `λ`.

For isotropic layers, returns `(n, n, n)` where `n` is the scalar refractive index.
For anisotropic layers, returns `(nx, ny, nz)` for each principal axis.
"""
function get_refractive_indices(layer::Layer, λ)
    if isanisotropic(layer)
        nx, ny, nz = layer.dispersion[1:3]
        return (nx(λ), ny(λ), nz(λ))
    else
        n = layer.dispersion(λ)
        return (n, n, n)
    end
end

"""
    get_euler_angles(layer::Layer)

Return the Euler angles `(φ, θ, ψ)` for crystal rotation.

For isotropic layers or unrotated anisotropic layers, returns `(0.0, 0.0, 0.0)`.
For rotated anisotropic layers, returns the stored angles.
"""
function get_euler_angles(layer::Layer)
    if isanisotropic(layer) && length(layer.dispersion) == 6
        return (layer.dispersion[4], layer.dispersion[5], layer.dispersion[6])
    else
        return (0.0, 0.0, 0.0)
    end
end

"""
    isrotated(layer::Layer)

Return `true` if the layer has non-zero Euler angles for crystal rotation.
"""
function isrotated(layer::Layer)
    φ, θ, ψ = get_euler_angles(layer)
    return φ != 0.0 || θ != 0.0 || ψ != 0.0
end

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
dielectric_constant(n::Real) = Complex(n^2)
