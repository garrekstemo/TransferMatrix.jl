struct Poynting
    out_p::SVector{3, Float64}
    in_p::SVector{3, Float64}
    out_s::SVector{3, Float64}
    in_s::SVector{3, Float64}
    refl_p::SVector{3, Float64}
    refl_s::SVector{3, Float64}

    function Poynting(out_p::T, in_p::T, out_s::T, in_s::T, refl_p::T, refl_s::T) where {T<:SVector{3, Float64}}
        new(out_p, in_p, out_s, in_s, refl_p, refl_s)
    end
end

"""
    ElectricField

Spatial electric-field profile through a layered structure, mirrored by
[`MagneticField`](@ref).

- `z`: position coordinates
- `p`: `(Ex, Ey, Ez)` for p-polarized incidence
- `s`: `(Ex, Ey, Ez)` for s-polarized incidence
- `boundaries`: z-positions of interfaces
"""
struct ElectricField{Z<:AbstractVector{Float64}}
    z::Z
    p::Matrix{ComplexF64}
    s::Matrix{ComplexF64}
    boundaries::Vector{Float64}

    function ElectricField(z::Z, p, s, boundaries) where {Z<:AbstractVector{Float64}}
        size(p, 2) == length(z) || throw(ArgumentError("p field columns must match z length"))
        size(s, 2) == length(z) || throw(ArgumentError("s field columns must match z length"))
        new{Z}(z, p, s, boundaries)
    end
end

"""
    MagneticField

Spatial magnetic-field profile through a layered structure, mirroring
[`ElectricField`](@ref). Fields are in impedance-normalized units `H̃ = Z₀ H_SI`
so `|E| ~ |H̃|` for a plane wave and E/H can be overlaid directly.

- `z`: position coordinates (same grid as `efield`)
- `p`: `(Hx, Hy, Hz)` for p-polarized incidence
- `s`: `(Hx, Hy, Hz)` for s-polarized incidence
- `boundaries`: z-positions of interfaces
"""
struct MagneticField{Z<:AbstractVector{Float64}}
    z::Z
    p::Matrix{ComplexF64}
    s::Matrix{ComplexF64}
    boundaries::Vector{Float64}

    function MagneticField(z::Z, p, s, boundaries) where {Z<:AbstractVector{Float64}}
        size(p, 2) == length(z) || throw(ArgumentError("p field columns must match z length"))
        size(s, 2) == length(z) || throw(ArgumentError("s field columns must match z length"))
        new{Z}(z, p, s, boundaries)
    end
end

"""
    TransferResult{T}

Container for reflectance and transmittance results from transfer matrix calculations.

# Fields
- `Tpp::T`: p-polarized transmittance (p-in → p-out)
- `Tss::T`: s-polarized transmittance (s-in → s-out)
- `Tps::T`: cross-polarized transmittance (p-in → s-out)
- `Tsp::T`: cross-polarized transmittance (s-in → p-out)
- `Rpp::T`: p-polarized reflectance (p-in → p-out)
- `Rss::T`: s-polarized reflectance (s-in → s-out)
- `Rps::T`: cross-polarized reflectance (p-in → s-out)
- `Rsp::T`: cross-polarized reflectance (s-in → p-out)

For single-wavelength calculations via `transfer()`, `T` is `Float64`.
For sweep calculations via `sweep_angle()` or `sweep_thickness()`, `T` is `Matrix{Float64}`.

!!! note "Cross-polarization terms"
    The cross-polarization terms (`Tps`, `Tsp`, `Rps`, `Rsp`) are zero for isotropic media
    and become non-zero for anisotropic (birefringent) materials.

# Examples
```julia
# Single wavelength - returns TransferResult{Float64}
result = transfer(1.0, layers)
result.Tpp  # Float64

# Sweep - returns TransferResult{Matrix{Float64}}
result = sweep_angle(λs, θs, layers)
result.Tpp  # Matrix{Float64}

# Destructuring works
(; Tpp, Rpp) = transfer(1.0, layers)
```
"""
struct TransferResult{T}
    Tpp::T
    Tss::T
    Tps::T
    Tsp::T
    Rpp::T
    Rss::T
    Rps::T
    Rsp::T
end

"""
    CircularTransferResult{T}

Circular-polarization (RCP/LCP) basis transmittance/reflectance — the analogue of
[`TransferResult`](@ref). Fields are right/left co- and cross-handedness terms:
`Trr` is `T_{R←R}`, `Trl` is `T_{R←L}`, etc.

As with [`TransferResult`](@ref), `T` is `Float64` for a single
`transfer(λ, layers; basis=:circular)` call and `Matrix{Float64}` for the sweep
functions.

# Convention
R/L are viewer-facing helicities defined in the fixed lab frame relative to `+z`,
under this package's `exp(-iωt)` time convention (basis states
`e_R = (x̂ - iŷ)/√2`, `e_L = (x̂ + iŷ)/√2`). Reflection flips the helicity label,
so at normal incidence on an isotropic interface the diagonal `Rrr`/`Rll` vanish
and reflection appears entirely in the off-diagonal `Rrl`/`Rlr`; at oblique
incidence the diagonal terms are small but nonzero. Unequal `Rrl ≠ Rlr` signals
genuine optical activity / magneto-optic coupling, not a bug.

See [`transfer`](@ref) for energy-ratio caveats.
"""
struct CircularTransferResult{T}
    Trr::T
    Tll::T
    Trl::T
    Tlr::T
    Rrr::T
    Rll::T
    Rrl::T
    Rlr::T
end
