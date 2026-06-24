"""
    Layer(material, thickness)
    Layer(nx, ny, nz, thickness)

Construct a single layer for transfer matrix calculations.

# Isotropic Layer (single refractive index)
```julia
Layer(material, thickness)
```
- `material`: a dispersion function `λ -> n(λ)`, or a `RefractiveMaterial` from
  RefractiveIndex.jl (the material form requires `using RefractiveIndex`, which
  activates the `RefractiveIndexExt` extension)
- `thickness`: Layer thickness in the same units as wavelength (typically μm)

# Anisotropic Layer (biaxial: three refractive indices)
```julia
Layer(nx, ny, nz, thickness)
```
- `nx`, `ny`, `nz`: Dispersion functions `λ -> n(λ)` for each principal axis
- `thickness`: Layer thickness in μm

For uniaxial materials, set two axes equal (e.g., `nx = ny` for optic axis along z).

# Units Convention
Lengths default to **micrometers (μm)** as bare numbers (`Layer(n, 0.1)` = 0.1 μm).
With `using Unitful`, thickness may carry units and is normalized to μm:
`Layer(n, 100u"nm")`. RefractiveIndex.jl dispersion functions are parameterized
in μm, so μm remains the internal unit.

`Layer` is parametric as `Layer{F,Mf,T}` where `F` is the dispersion function type, `Mf` is the
permeability type (`Nothing` or a callable `λ -> SMatrix{3,3,ComplexF64}`), and `T` is the thickness type.
For anisotropic layers, `F` is a `Tuple` of three dispersion functions.

# Examples
```julia
# Isotropic: Using a RefractiveIndex.jl material (requires `using RefractiveIndex`)
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
struct Layer{F,Mf,T<:Real}
    dispersion::F
    mu::Mf
    thickness::T

    function Layer(dispersion::F, mu::Mf, thickness::T) where {F,Mf,T<:Real}
        thickness ≥ 0 || throw(DomainError("Layer thickness must be non-negative"))
        new{F,Mf,T}(dispersion, mu, thickness)
    end
end

# Normalize the `mu=` input to `nothing` or a callable λ -> SMatrix{3,3,ComplexF64}.
_normalize_mu(::Nothing) = nothing
_normalize_mu(m::Number) = let M = SMatrix{3,3,ComplexF64}(m*I); λ -> M end
_normalize_mu(M::AbstractMatrix) = let Ms = SMatrix{3,3,ComplexF64}(M); λ -> Ms end
_normalize_mu(f) = λ -> SMatrix{3,3,ComplexF64}(f(λ))

# Isotropic constructor: material is any callable λ -> n(λ).
# The Layer(::RefractiveMaterial, thickness) constructor lives in the RefractiveIndexExt
# package extension (load RefractiveIndex to enable it).
Layer(material, thickness::Real; mu=nothing) =
    Layer(material, _normalize_mu(mu), thickness)

# Isotropic constructor from tabulated (λ, n, k) data.
Layer(λs::AbstractVector, dispersion::AbstractVector, extinction::AbstractVector, thickness::Real; mu=nothing) =
    Layer(refractive_index(λs, dispersion, extinction), thickness; mu=mu)

"""
    Base.broadcastable(layer::Layer)

Treat a single `Layer` as a scalar under broadcasting so helpers can be broadcast
over one layer across many wavelengths without an explicit `Ref`, e.g.
`get_refractive_indices.(layer, λs)`. Without this, broadcasting a bare `Layer`
hits Julia's `collect`-based fallback and throws a confusing
`MethodError: no method matching length(::Layer)`.

A `Vector{Layer}` stack still broadcasts element-wise, which is the desired
behavior for per-layer operations. For threaded wavelength/angle/thickness sweeps,
prefer the dedicated `sweep_angle` / `sweep_thickness` API.
"""
Base.broadcastable(layer::Layer) = Ref(layer)

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
function Layer(nx::F1, ny::F2, nz::F3, thickness::T; euler::NTuple{3,Real}=(0.0, 0.0, 0.0), mu=nothing) where {F1,F2,F3,T<:Real}
    thickness ≥ 0 || throw(DomainError("Layer thickness must be non-negative"))
    φ, θ, ψ = Float64.(euler)
    Layer((nx, ny, nz, φ, θ, ψ), _normalize_mu(mu), thickness)
end

"""
    isanisotropic(layer::Layer)

Return `true` if the layer has anisotropic optical properties (different refractive
indices along principal axes), `false` for isotropic layers.
"""
isanisotropic(layer::Layer) = layer.dispersion isa Tuple

"""
    ismagnetic(layer::Layer)

Return `true` if the layer has a non-trivial magnetic permeability (mu ≠ nothing).
"""
ismagnetic(layer::Layer) = layer.mu !== nothing

"""
    get_permeability(layer::Layer, λ)

Return the 3×3 permeability tensor `SMatrix{3,3,ComplexF64}` for the layer at
wavelength `λ`, or `nothing` if the layer is non-magnetic.
"""
get_permeability(layer::Layer, λ) = layer.mu === nothing ? nothing : layer.mu(λ)

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
    refractive_index(λs, ns, ks)

Build a complex dispersion function `λ -> n(λ) + i·k(λ)` by linearly interpolating
tabulated real (`ns`) and imaginary (`ks`) refractive-index data over wavelengths `λs`.

A `refractive_index(material::RefractiveMaterial)` method that derives the dispersion
function from a RefractiveIndex.jl material is provided by the `RefractiveIndexExt`
package extension — load `RefractiveIndex` to enable it (along with the
`Layer(::RefractiveMaterial, d)` and `Sheet(::RefractiveMaterial, d)` constructors).
"""
function refractive_index(λs::AbstractVector, ns::AbstractVector, ks::AbstractVector)
    n_real = LinearInterpolation(ns, λs)
    n_imag = LinearInterpolation(ks, λs)
    return λ -> begin
        n_real(λ) + im * n_imag(λ)
    end
end

"""
    find_bounds(layers)

Return cumulative thickness positions measured from the start of layer 1, along with
the total thickness.

The returned vector has one entry per layer, where entry `i` is the sum of thicknesses
from layer 1 through layer `i`. These are **not** interface positions relative to z = 0;
callers that need z-coordinates (e.g., `efield`) must subtract the first layer thickness.
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
