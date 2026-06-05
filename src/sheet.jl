"""
    Sheet(σ)
    Sheet(; xx, yy, xy=0, yx=0)
    Sheet(material, d)
    Sheet(nx, ny, d)

A zero-thickness 2D conductive sheet (e.g. a TMDC monolayer) for the transfer
matrix. The internal representation is a callable `λ -> SMatrix{2,2,ComplexF64}`
returning the **SI sheet conductivity tensor** `[σxx σxy; σyx σyy]` in Siemens.

- `Sheet(σ)` — scalar isotropic conductivity (`σ::Number` or `σ(λ)->Number`),
  stored as `diag(σ, σ)`.
- `Sheet(; xx, yy, xy=0, yx=0)` — anisotropic tensor; each entry a `Number` or
  `λ->Number`.
- `Sheet(material, d)` — convert a refractive index to a sheet via
  `σ = -i ω ε₀ d (n²-1)` (isotropic). `material` is a `RefractiveMaterial` or
  `n(λ)`; `d` is the effective thickness in μm.
- `Sheet(nx, ny, d)` — in-plane anisotropic index conversion (`σxx` from `nx`,
  `σyy` from `ny`).

The factor `Z₀ = √(μ₀/ε₀)` is applied exactly once, later, in [`sheet_matrix`](@ref).
"""
struct Sheet{F}
    conductivity::F
    Sheet{F}(c::F) where {F} = new{F}(c)
end

# Internal: store a raw λ -> SMatrix closure without scalar-wrapping.
_rawsheet(f) = Sheet{typeof(f)}(f)

# Internal: evaluate a component that may be a constant or a function of λ.
_cval(x::Number, λ) = ComplexF64(x)
_cval(f, λ) = ComplexF64(f(λ))

# Internal: turn an index spec into an n(λ) function.
_index_fn(m::RefractiveMaterial) = refractive_index(m)
_index_fn(x::Number) = _ -> x
_index_fn(f::Function) = f

# Internal: SI scalar sheet conductivity from an index function and thickness d (μm).
_sigma_from_index(nfun, d) = λ -> -im * (2π * c_0 / λ) * ε_0 * d * (ComplexF64(nfun(λ))^2 - 1)

# Internal: wrap a scalar conductivity (constant or function) as a diagonal sheet.
_diagonal_sheet(σ) = _rawsheet(λ -> (s = _cval(σ, λ); SMatrix{2,2,ComplexF64}(s, 0, 0, s)))

# Scalar isotropic conductivity (Number or function returning a Number).
Sheet(σ::Number) = _diagonal_sheet(σ)
Sheet(σ::Function) = _diagonal_sheet(σ)

# Anisotropic tensor by component. SMatrix is column-major: (xx, yx, xy, yy).
function Sheet(; xx, yy, xy = 0, yx = 0)
    _rawsheet(λ -> SMatrix{2,2,ComplexF64}(_cval(xx, λ), _cval(yx, λ), _cval(xy, λ), _cval(yy, λ)))
end

# Index -> conductivity (isotropic).
Sheet(material::RefractiveMaterial, d::Real) = _diagonal_sheet(_sigma_from_index(_index_fn(material), d))
Sheet(n::Function, d::Real) = _diagonal_sheet(_sigma_from_index(n, d))

# Index -> conductivity (in-plane anisotropic).
function Sheet(nx, ny, d::Real)
    σx = _sigma_from_index(_index_fn(nx), d)
    σy = _sigma_from_index(_index_fn(ny), d)
    _rawsheet(λ -> SMatrix{2,2,ComplexF64}(σx(λ), 0, 0, σy(λ)))
end
