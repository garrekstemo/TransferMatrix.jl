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
permeability type (`Nothing` or a callable `λ -> SMatrix{3,3,<:Complex}`), and `T` is the thickness type.
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

# Magnetic / Permeability (`mu=`)

All constructors accept an optional `mu=` keyword that sets the layer's magnetic
permeability. Accepted forms: `nothing` (default, use the global `μ=` fallback),
a scalar `Number` (isotropic μ), a constant `3×3 AbstractMatrix`, or a callable
`λ -> 3×3 matrix` for wavelength-dependent tensors.

```julia
# Gyromagnetic (Polder) ferrite film, bias ∥ z
film = Layer(λ -> 1.0, 0.4; mu = gyrotropic_tensor(2.0, 0.6))

# Uniaxial-μ slab
slab = Layer(λ -> 1.5, 0.3; mu = [2 0 0; 0 2 0; 0 0 3])
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
_normalize_mu(m::Number) = let M = SMatrix{3,3,complex(float(typeof(m)))}(m*I); λ -> M end
_normalize_mu(M::AbstractMatrix) = let Ms = _complex_smatrix3(M); λ -> Ms end
_normalize_mu(f) = λ -> _complex_smatrix3(f(λ))

# 3×3 complex SMatrix with the element type of `M` (ComplexF64 for Float64 input, a
# complex dual number under ForwardDiff), so a fitted μ tensor stays differentiable.
_complex_smatrix3(M) = SMatrix{3,3,complex(float(eltype(M)))}(M)

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
    φ, θ, ψ = float.(promote(euler...))
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

Return the 3×3 complex permeability tensor `SMatrix{3,3}` for the layer at
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

# Immutable linear interpolation over sorted knots. Stateless by design: dispersion
# closures are evaluated concurrently by the threaded sweeps (`sweep_angle` /
# `sweep_thickness` share one `layers` vector across threads), so evaluation must
# not mutate the interpolant. DataInterpolations' LinearInterpolation writes a
# search-hint `Base.RefValue` on every call, which was a data race here.
struct TabulatedInterpolation{V<:AbstractVector, W<:AbstractVector}
    knots::V
    values::W

    function TabulatedInterpolation(knots::V, values::W) where {V<:AbstractVector, W<:AbstractVector}
        length(knots) == length(values) ||
            throw(ArgumentError("knot and value vectors must have equal length"))
        length(knots) >= 2 ||
            throw(ArgumentError("interpolation requires at least 2 data points"))
        # Strictly increasing: a duplicated knot creates a zero-width segment,
        # which evaluates to NaN if λ lands on it (e.g. a duplicated endpoint).
        issorted(knots, lt = <=) ||
            throw(ArgumentError("wavelength knots must be strictly increasing (no duplicates)"))
        new{V,W}(knots, values)
    end
end

function (itp::TabulatedInterpolation)(λ)
    knots, values = itp.knots, itp.values
    first(knots) <= λ <= last(knots) ||
        throw(DomainError(λ, "wavelength outside tabulated range [$(first(knots)), $(last(knots))]"))
    i = min(searchsortedlast(knots, λ), lastindex(knots) - 1)
    t = (λ - knots[i]) / (knots[i+1] - knots[i])
    return values[i] + t * (values[i+1] - values[i])
end

"""
    refractive_index(λs, ns, ks)

Build a complex dispersion function `λ -> n(λ) + i·k(λ)` by linearly interpolating
tabulated real (`ns`) and imaginary (`ks`) refractive-index data over wavelengths `λs`
(sorted, increasing). Evaluation outside the tabulated range throws a `DomainError`.

A `refractive_index(material::RefractiveMaterial)` method that derives the dispersion
function from a RefractiveIndex.jl material is provided by the `RefractiveIndexExt`
package extension — load `RefractiveIndex` to enable it (along with the
`Layer(::RefractiveMaterial, d)` and `Sheet(::RefractiveMaterial, d)` constructors).
"""
function refractive_index(λs::AbstractVector, ns::AbstractVector, ks::AbstractVector)
    n_real = TabulatedInterpolation(λs, ns)
    n_imag = TabulatedInterpolation(λs, ks)
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
    Td = float(mapreduce(l -> typeof(l.thickness), promote_type, layers; init=Float64))
    interface_positions = Vector{Td}(undef, n_layers)

    total_thickness = zero(Td)
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


# ---------------------------------------------------------------------------
# Number type of a calculation.
#
# Every propagation core works in one complex element type, `Complex{R}`, where R
# is the promotion of the real types of λ, θ, the global μ, each layer's thickness,
# refractive indices, Euler angles and permeability, and each sheet's conductivity.
# For ordinary inputs R = Float64 (so the cores keep their ComplexF64 code paths);
# under ForwardDiff R is a Dual, which is what lets gradients pass through.
#
# Refractive indices are the one input whose type is only known by evaluating the
# dispersion function, which can be the most expensive step of a call. So `transfer`
# first runs with the type of everything else (`_stack_realtype_static`), and each
# function barrier that evaluates a dispersion checks that the value fits
# (`_fits`). Only if one does not (a dual-valued index while everything else is
# Float64, i.e. when fitting indices) is the point recomputed with the full type
# (`_stack_realtype`), so ordinary calls evaluate each dispersion once.

_realtype(::Type{T}) where {T<:Number} = real(T)
_realtype(::Type{T}) where {T<:AbstractArray} = real(eltype(T))
_realtype(::Type) = Bool

_mu_realtype(::Nothing, λ) = Bool
_mu_realtype(f, λ) = _realtype(typeof(f(λ)))

_euler_realtype(t::Tuple) = length(t) == 6 ? promote_type(map(typeof, t[4:6])...) : Bool
_euler_realtype(_) = Bool

_layer_static_realtype(layer::Layer{F,Mf,Td}, λ) where {F,Mf,Td} =
    promote_type(Td, _euler_realtype(layer.dispersion), _mu_realtype(layer.mu, λ))

_dispersion_realtype(f, λ) = _realtype(typeof(f(λ)))
_dispersion_realtype(t::Tuple, λ) =
    promote_type(_dispersion_realtype(t[1], λ), _dispersion_realtype(t[2], λ), _dispersion_realtype(t[3], λ))

_layer_realtype(layer::Layer, λ) =
    promote_type(_layer_static_realtype(layer, λ), _dispersion_realtype(layer.dispersion, λ))

function _accumulate_realtype(f::F, λ, θ, μ, layers, sheets) where {F}
    R = float(promote_type(_realtype(typeof(λ)), _realtype(typeof(θ)), _realtype(typeof(μ))))
    for layer in layers
        Rl = f(layer, λ)::DataType
        Rl === R || (R = promote_type(R, Rl))
    end
    if sheets !== nothing
        for sheet in values(sheets)
            Rs = _realtype(typeof(sheet.conductivity(λ)))
            Rs === R || (R = promote_type(R, Rs))
        end
    end
    return R
end

# Real type of everything except the refractive indices (no dispersion evaluated).
_stack_realtype_static(λ, θ, μ, layers, sheets=nothing) =
    _accumulate_realtype(_layer_static_realtype, λ, θ, μ, layers, sheets)

# Full real type of the calculation (evaluates every dispersion once).
_stack_realtype(λ, θ, μ, layers, sheets=nothing) =
    _accumulate_realtype(_layer_realtype, λ, θ, μ, layers, sheets)

"""
    _stack_ctype(λ, θ, μ, layers, sheets=nothing)

Complex element type used by the propagation cores for this point: `ComplexF64` for
ordinary inputs, `Complex{<:ForwardDiff.Dual}` when any input carries dual numbers.
"""
_stack_ctype(λ, θ, μ, layers, sheets=nothing) = Complex{_stack_realtype(λ, θ, μ, layers, sheets)}

# Whether a value computed from a dispersion function can be stored in element type C.
@inline function _fits(::Type{C}, x) where {C}
    Rx = _realtype(typeof(x))
    return Rx === real(C) || promote_type(real(C), Rx) === real(C)
end
@inline _fits(::Type{C}, x, y, z) where {C} = _fits(C, x) && _fits(C, y) && _fits(C, z)

# Returned by a propagation core run with a too-narrow element type (see above).
struct _Widen end
const _WIDEN = _Widen()
