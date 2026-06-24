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


# Normalize accepted sheet inputs (Dict or iterable of `i => sheet` pairs) to Dict{Int,Sheet}.
_sheets_dict(s::Dict{Int,Sheet}) = s
_sheets_dict(s) = Dict{Int,Sheet}(Int(k) => v for (k, v) in s)

# Unconditional structural validation: keys must index an interior interface.
function _validate_sheet_indices(sd, N)
    sd === nothing && return nothing
    for i in keys(sd)
        (1 ≤ i ≤ N - 1) || throw(ArgumentError("sheet index $i out of range; must be 1 ≤ i ≤ $(N - 1)"))
    end
    return nothing
end


"""
    poynting(Ψ, a)

Calculates the Poynting vector for the structure
from Ψ and matrix ``a``.

From Berreman, 1972, Ψ is the column matrix:

```math
\\Psi = 
    \\begin{pmatrix}
        Ex \\\\\
        Hy \\\\\
        Ey \\\\\
       -Hx
    \\end{pmatrix}
```

for a right-handed Cartesian coordinate system with
the z-axis along the normal to the multilayer structure.

Berreman, 1972,
https://doi.org/10.1364/JOSA.62.000502
"""
function poynting(Ψ, a)

    S = @MMatrix zeros(ComplexF64, 3, 4)

    for m in 1:4
        Ex =  Ψ[1, m]
        Ey =  Ψ[3, m]
        Hx = -Ψ[4, m]
        Hy =  Ψ[2, m]

        Ez = a[3,1] * Ex + a[3,2] * Ey + a[3,4] * Hx + a[3,5] * Hy
        Hz = a[6,1] * Ex + a[6,2] * Ey + a[6,4] * Hx + a[6,5] * Hy
        
        S[1, m] = Ey * Hz - Ez * Hy
        S[2, m] = Ez * Hx - Ex * Hz
        S[3, m] = Ex * Hy - Ey * Hx
    end
    return SMatrix(S)
end


"""
    poynting(ξ, q_in, q_out, γ_in, γ_out, t_coefs, r_coefs[, μ_in, μ_out])

Calculate the Poynting vector from wavevectors ``q``,
components of the electric field γ, and transmission
and reflection coefficients.

Transmitted Poynting vectors use substrate wavevectors (`q_out`), while
reflected Poynting vectors use incident-medium wavevectors (`k_in[3,:]`,
`k_in[4,:]`), since reflected waves propagate in the incident medium.

The optional `μ_in` and `μ_out` arguments are the 3×3 permeability tensors of
the ambient (incident) and substrate media respectively (default: identity, i.e.
non-magnetic). They are used to form ``H = μ^{-1}(k̄ \\times E)`` when computing
the Poynting flux, so that energy conservation holds for a magnetic substrate.
A magnetic ambient is only correctly handled at normal incidence; at oblique
incidence the conserved in-plane wavevector ξ is still computed from the
ambient permittivity alone (see issue #71), so results may be inaccurate.

!!! note "Transmittance vs reflectance"
    This function computes Poynting vectors for both transmitted and reflected
    waves, but only the **transmitted** Poynting vectors are used for the final
    output. Reflectance is computed as ``R = |r|^2`` from the transfer matrix
    coefficients — see [`transfer`](@ref) for the rationale.
"""
function poynting(ξ, q_in, q_out, γ_in, γ_out, t_coefs, r_coefs,
                  μ_in::AbstractMatrix  = SMatrix{3,3,ComplexF64}(I),
                  μ_out::AbstractMatrix = SMatrix{3,3,ComplexF64}(I))

    μin_inv  = inv(SMatrix{3,3,ComplexF64}(μ_in))
    μout_inv = inv(SMatrix{3,3,ComplexF64}(μ_out))

    # create the wavevector in the first layer
    k_in = @MMatrix zeros(ComplexF64, 4, 3)
    k_in[:, 1] .= ξ

    for (i, q_i) in enumerate(q_in)
        k_in[i, 3] = q_i
    end

    k_in ./= c_0
    k_in = SMatrix(k_in)

    E_forward_in_p =  γ_in[1, :]  # p-polarized incident electric field
    E_forward_in_s =  γ_in[2, :]  # s-polarized incident electric field
    # E_backward_in_p = γ_in[3, :]
    # E_backward_in_s = γ_in[4, :]

    # Each transmitted/reflected field is a superposition of the two
    # substrate (resp. incident) eigenmodes, which carry the substrate mode-1
    # field γ_out[1] and mode-2 field γ_out[2].
    E_out_p1 = t_coefs[1] * γ_out[1, :]
    E_out_p2 = t_coefs[2] * γ_out[2, :]
    E_out_s1 = t_coefs[3] * γ_out[1, :]
    E_out_s2 = t_coefs[4] * γ_out[2, :]

    E_ref_p = r_coefs[1] * γ_in[3, :] + r_coefs[2] * γ_in[4, :]
    E_ref_s = r_coefs[3] * γ_in[3, :] + r_coefs[4] * γ_in[4, :]

    S_in_p = real(0.5 * E_forward_in_p × conj(μin_inv * (k_in[1, :] × E_forward_in_p)))
    S_in_s = real(0.5 * E_forward_in_s × conj(μin_inv * (k_in[2, :] × E_forward_in_s)))

    k_out = @MMatrix zeros(ComplexF64, 4, 3)
    k_out[:, 1] .= ξ

    for (i, q_i) in enumerate(q_out)
        k_out[i, 3] = q_i
    end

    k_out ./= c_0
    k_out = SMatrix(k_out)

    # Transmitted power is the sum of the two substrate eigenmodes' Poynting
    # vectors, EACH evaluated with its OWN wavevector (k_out modes 1, 2). The
    # two modes carry different wavevectors when the substrate is anisotropic,
    # and the cross-interference between distinct propagating modes transports
    # no net z-flux, so the per-mode sum is the physical transmitted flux.
    # (For an isotropic substrate k_out[1]=k_out[2] and this reduces to the
    # single-wavevector form, since the p/s cross term carries no z-power.)
    S_out_p = real(0.5 * E_out_p1 × conj(μout_inv * (k_out[1, :] × E_out_p1))) +
              real(0.5 * E_out_p2 × conj(μout_inv * (k_out[2, :] × E_out_p2)))
    S_out_s = real(0.5 * E_out_s1 × conj(μout_inv * (k_out[1, :] × E_out_s1))) +
              real(0.5 * E_out_s2 × conj(μout_inv * (k_out[2, :] × E_out_s2)))

    # Reflected waves: use incident-medium wavevectors (k_in modes 3,4),
    # not substrate wavevectors, because reflected light propagates backward
    # through the incident medium with its own q-values.
    S_refl_p = real(0.5 * E_ref_p × conj(μin_inv * (k_in[3, :] × E_ref_p)))
    S_refl_s = real(0.5 * E_ref_s × conj(μin_inv * (k_in[4, :] × E_ref_s)))

    return Poynting(S_out_p, S_in_p, S_out_s, S_in_s, S_refl_p, S_refl_s)
end


"""
    evaluate_birefringence(Ψ, S, t_modes, r_modes)

For the four modes (two transmitting and two reflecting), the ratio

```math
\\begin{aligned}
    C &= |E_x|^2 / (|E_x|^2 + |E_y|^2) \\\\\
      &= |Ψ_1|^2 / (|Ψ_1|^2 + |Ψ_3|^2)
\\end{aligned}
```

is evaluated. Recall that the values for the electric field are contained
in the eigenvector matrix, Ψ.

If the layer material is birefringent, there will be anisotropy in the
dielectric tensor. If this is the case, the x and y components of the 
Poynting vector needs to be analyzed (eqn 15 in Passler et al., 2017):

```math
C = |S_x|^2 / (|S_x|^2 + |S_y|^2)
```

If there is no birefringence, then the electric field is analyzed.
This analysis follows

Li et al., 1988, https://doi.org/10.1364/AO.27.001334

and the use of the Poynting vector is from

Passler et al., 2017, https://doi.org/10.1364/JOSAB.34.002128
Passler et al., 2019, https://doi.org/10.1364/JOSAB.36.003246
"""
# Sort one mode pair (a transmitted or reflected pair) so the p-like mode is
# first and the s-like mode second, matching the fixed polarization references
# in `calculate_γ`.
#
# The Poynting ratio C = |Sx|²/(|Sx|²+|Sy|²) distinguishes p from s only when
# the modes actually differ along x vs y. For an axis-aligned (diagonal-ε)
# crystal BOTH eigenmodes have Sy = 0, so C ≈ 1 for both and the ratio is
# degenerate — it cannot tell p from s. In that case fall back to the
# electric-field ratio |Ex|²/(|Ex|²+|Ey|²), which separates the p-mode
# (Ex ≠ 0) from the s-mode (Ey ≠ 0). `isapprox(NaN, NaN) = false`, so a 0/0
# Poynting ratio (both Sx and Sy zero) is also routed to the E-field fallback
# via the `!isfinite` checks.
function sort_polarization_pair!(modes, Ψ, S)
    Cp = abs_ratio(S[1, modes[1]], S[2, modes[1]])
    Cs = abs_ratio(S[1, modes[2]], S[2, modes[2]])
    if isapprox(Cp, Cs) || !isfinite(Cp) || !isfinite(Cs)
        Cp = abs_ratio(Ψ[1, modes[1]], Ψ[3, modes[1]])
        Cs = abs_ratio(Ψ[1, modes[2]], Ψ[3, modes[2]])
    end
    if Cs > Cp
        reverse!(modes)
    end
    return modes
end

function evaluate_birefringence(Ψ, S, t_modes, r_modes)

    # Sort each pair (transmitted, reflected) so the p-like mode is first
    # (slot 1 / slot 3) and the s-like mode second (slot 2 / slot 4). This
    # ordering is assumed by the fixed polarization references in
    # `calculate_γ` (γ[1,1]=1, γ[2,2]=1, γ[3,1]=-1, γ[4,2]=1); getting it
    # wrong assigns each mode the wrong polarization, which both makes the
    # cross-polarization denominators in `calculate_γ` vanish (0/0 → NaN) and
    # corrupts the r/t coefficients (energy is not conserved).
    sort_polarization_pair!(t_modes, Ψ, S)
    sort_polarization_pair!(r_modes, Ψ, S)

    return t_modes, r_modes
end


"""
Ratio of the absolution squares of two components
used to evaluate if a material is birefringent.
"""
abs_ratio(a, b) = abs2(a) / (abs2(a) + abs2(b))


# Lightweight path: only computes Γ and S without allocating per-layer
# Ds, Ps, γs vectors. Used by `transfer` in tight spectral loops.
function _propagate_core(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

    N = length(layers)
    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ε_0in = dielectric_constant(nx_in)
    ξ = √(ε_0in) * sin(θ)

    Λ_1324 = @SMatrix [1 0 0 0;
                       0 0 1 0;
                       0 1 0 0;
                       0 0 0 1]

    no_sheets = sheets === nothing || isempty(sheets)

    D_prev, _, γ_first, q_first = layer_matrices(layers[1], λ, ξ, μ)
    γ_last = γ_first
    q_last = q_first

    Γ = SMatrix{4,4,ComplexF64}(I)
    for i in 2:N
        layer = layers[i]
        D_cur, P_cur, γ_cur, q_cur = layer_matrices(layer, λ, ξ, μ)
        if no_sheets || !haskey(sheets, i - 1)
            L = D_prev \ D_cur                                  # interface (i-1, i)
        else
            L = D_prev \ (sheet_matrix(sheets[i - 1], λ) * D_cur)
        end
        Γ *= L                                                  # first ⇒ D₀⁻¹D₂ ; last ⇒ D_{N-1}⁻¹D_f
        if i < N
            Γ *= P_cur(layer.thickness)                         # propagate interior layer i
        end
        D_prev = D_cur
        if i == N
            γ_last = γ_cur
            q_last = q_cur
        end
    end

    Γ = (Λ_1324 \ Γ) * Λ_1324
    r, R, t, T = calculate_tr(Γ)
    μ_in_mat  = ismagnetic(layers[1])   ? get_permeability(layers[1],   λ) : SMatrix{3,3,ComplexF64}(μ*I)
    μ_out_mat = ismagnetic(layers[end]) ? get_permeability(layers[end], λ) : SMatrix{3,3,ComplexF64}(μ*I)
    S = poynting(ξ, q_first, q_last, γ_first, γ_last, t, r, μ_in_mat, μ_out_mat)

    return Γ, S
end


"""
    _propagate_core_exp(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Matrix-exponential propagation core. Interior layers propagate via
[`layer_transfer_exp`](@ref) (no eigenmode sorting); the semi-infinite ambient and
substrate keep the eigenmode treatment in [`layer_matrices`](@ref), which is needed
for the r/t coefficients and the Poynting transmittance. Conductive sheets are
injected at their respective interfaces, mirroring [`_propagate_core`](@ref).
Returns `(Γ, S)` like [`_propagate_core`](@ref).
"""
function _propagate_core_exp(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

    N = length(layers)
    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ε_0in = dielectric_constant(nx_in)
    ξ = √(ε_0in) * sin(θ)
    ω = 2π * c_0 / λ

    D_1, _, γ_first, q_first = layer_matrices(layers[1], λ, ξ, μ)
    D_N, _, γ_last,  q_last  = layer_matrices(layers[N], λ, ξ, μ)

    # Interior product in the dynamical-matrix basis. Tangential fields are
    # continuous across plain interfaces, so interior layers chain directly; a
    # conductive sheet at interface (i, i+1) is injected as sheet_matrix.
    no_sheets = sheets === nothing || isempty(sheets)
    core = SMatrix{4,4,ComplexF64}(I)
    if !no_sheets && haskey(sheets, 1)
        core = core * sheet_matrix(sheets[1], λ)
    end
    for i in 2:N-1
        core = core * layer_transfer_exp(layers[i], λ, ξ, ω, μ)
        if !no_sheets && haskey(sheets, i)
            core = core * sheet_matrix(sheets[i], λ)
        end
    end
    core = core * D_N

    Γ = (_Λ1324 \ (D_1 \ core)) * _Λ1324
    r, R, t, T = calculate_tr(Γ)
    μ_in_mat  = ismagnetic(layers[1])   ? get_permeability(layers[1],   λ) : SMatrix{3,3,ComplexF64}(μ*I)
    μ_out_mat = ismagnetic(layers[end]) ? get_permeability(layers[end], λ) : SMatrix{3,3,ComplexF64}(μ*I)
    S = poynting(ξ, q_first, q_last, γ_first, γ_last, t, r, μ_in_mat, μ_out_mat)

    return Γ, S
end

# Dispatch between the eigenmode (:eig) and matrix-exponential (:exp) cores.
_propagate(::Val{:exp}, λ, layers; kwargs...) = _propagate_core_exp(λ, layers; kwargs...)
_propagate(::Val{:eig}, λ, layers; kwargs...) = _propagate_core(λ, layers; kwargs...)
_propagate(::Val{M}, λ, layers; kwargs...) where {M} =
    throw(ArgumentError("method must be :exp or :eig, got :$(M)"))


"""
    _propagate_full(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Internal full transfer-matrix pass. Returns `(Γ, S, Ds, Ps, γs, qs)` — like
[`propagate`](@ref) but also returns the per-layer eigenvalue vectors `qs`,
needed for magnetic-field reconstruction. Supports conductive sheets (Task: sheets).
"""
function _propagate_full(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

    N = length(layers)
    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ε_0in = dielectric_constant(nx_in)
    ξ = √(ε_0in) * sin(θ)

    Λ_1324 = @SMatrix [1 0 0 0;
                       0 0 1 0;
                       0 1 0 0;
                       0 0 0 1]

    no_sheets = sheets === nothing || isempty(sheets)

    D_1, P_1, γ_1, q_1 = layer_matrices(layers[1], λ, ξ, μ)
    Ds = Vector{typeof(D_1)}(undef, N)
    Ps = Vector{typeof(P_1)}(undef, N)
    γs = Vector{typeof(γ_1)}(undef, N)
    qs = Vector{typeof(q_1)}(undef, N)
    Ds[1] = D_1; Ps[1] = P_1; γs[1] = γ_1; qs[1] = q_1

    Γ = SMatrix{4,4,ComplexF64}(I)
    for i in 2:N
        D_i, P_i, γ_i, q_i = layer_matrices(layers[i], λ, ξ, μ)
        Ds[i] = D_i; Ps[i] = P_i; γs[i] = γ_i; qs[i] = q_i
        if no_sheets || !haskey(sheets, i - 1)
            L = Ds[i - 1] \ D_i
        else
            L = Ds[i - 1] \ (sheet_matrix(sheets[i - 1], λ) * D_i)
        end
        Γ *= L
        if i < N
            Γ *= P_i(layers[i].thickness)
        end
    end

    Γ = (Λ_1324 \ Γ) * Λ_1324
    r, R, t, T = calculate_tr(Γ)
    μ_in_mat  = ismagnetic(layers[1])   ? get_permeability(layers[1],   λ) : SMatrix{3,3,ComplexF64}(μ*I)
    μ_out_mat = ismagnetic(layers[end]) ? get_permeability(layers[end], λ) : SMatrix{3,3,ComplexF64}(μ*I)
    S = poynting(ξ, q_1, qs[N], γ_1, γs[N], t, r, μ_in_mat, μ_out_mat)

    return Γ, S, Ds, Ps, γs, qs
end

"""
    propagate(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Calculate the transfer matrix and Poynting vector for the structure, plus the
per-layer `D`, `P`, and `γ` matrices used for field reconstruction. Returns the
5-tuple `(Γ, S, Ds, Ps, γs)`. See [`transfer`](@ref) for the public R/T API.
"""
propagate(λ, layers; θ=0.0, μ=1.0, sheets=nothing) =
    _propagate_full(λ, layers; θ=θ, μ=μ, sheets=sheets)[1:5]


# Linear (p,s) -> circular (L,R) change-of-basis matrices for the exp(-iωt) time
# convention. Columns of C are ordered (L, R): e_L = (x̂ + iŷ)/√2, e_R = (x̂ - iŷ)/√2,
# so [Ex; Ey] = C [E_L; E_R]. The sign of i is fixed by exp(-iωt) (do NOT flip it
# without also changing the package time convention).
const _C_CIRC = SMatrix{2,2,ComplexF64}(1, im, 1, -im) / sqrt(2)      # [1 1; im -im]
const _C_CIRC_INV = SMatrix{2,2,ComplexF64}(1, 1, -im, im) / sqrt(2)  # [1 -im; 1 im]

# Transmission-null tolerance: below this the transmitted Poynting flux is taken
# as zero, so circular transmittance is defined as 0 (avoids 0/0 at a perfect mirror).
const _CIRC_T_TOL = 1e-12

# Transform a 2×2 linear Jones matrix ([out,in] ordered p,s) into the circular
# (L,R) basis: J_circ = C⁻¹ J_lin C. Returns a 2×2 ([out,in] ordered L,R).
_jones_to_circular(J_lin) = _C_CIRC_INV * J_lin * _C_CIRC

# Assemble circular-basis R/T from the linear complex Jones coefficients
# r = (rpp,rps,rss,rsp), t = (tpp,tps,tsp,tss) and the Poynting diagonal
# transmittances Tpp, Tss. Reflectance is |r_circ|^2; transmittance is N·|t_circ|^2
# with the single polarization-independent Poynting scalar N (exact for an isotropic
# substrate; see `transfer` docstring for caveats).
function _circular_result(r, t, Tpp, Tss)
    r_lin = SMatrix{2,2,ComplexF64}(r[1], r[4], r[2], r[3])  # [rpp rps; rsp rss]
    t_lin = SMatrix{2,2,ComplexF64}(t[1], t[3], t[2], t[4])  # [tpp tps; tsp tss]
    r_c = _jones_to_circular(r_lin)
    t_c = _jones_to_circular(t_lin)
    Rc = abs2.(r_c)
    denom = abs2(t[1]) + abs2(t[4])  # |tpp|^2 + |tss|^2
    N = denom > _CIRC_T_TOL ? (Tpp + Tss) / denom : zero(denom)
    Tc = N .* abs2.(t_c)
    # circular matrices are [out,in] ordered (L,R): [1,1]=LL [1,2]=LR [2,1]=RL [2,2]=RR
    return CircularTransferResult(Tc[2, 2], Tc[1, 1], Tc[2, 1], Tc[1, 2],
                                  Rc[2, 2], Rc[1, 1], Rc[2, 1], Rc[1, 2])
end


"""
    calculate_tr(Γ)

Calculate reflectance and transmittance for the total stack.
This takes the matrix Γ* in Passler, et al., but for brevity we call it Γ in this function.

The original formalism is from:
Yeh, 1979,
https://doi.org/10.1364/JOSA.69.000742

but the ordering of reflection/transmission coefficients 
is modified in Passler, et al. 2017
https://doi.org/10.1364/JOSAB.34.002128
"""
function calculate_tr(Γ)

    d = Γ[1,1] * Γ[3,3] - Γ[1,3] * Γ[3,1]

    rpp = (Γ[2,1] * Γ[3,3] - Γ[2,3] * Γ[3,1]) / d
    rss = (Γ[1,1] * Γ[4,3] - Γ[4,1] * Γ[1,3]) / d
    rps = (Γ[4,1] * Γ[3,3] - Γ[4,3] * Γ[3,1]) / d
    rsp = (Γ[1,1] * Γ[2,3] - Γ[2,1] * Γ[1,3]) / d

    tpp =  Γ[3,3] / d
    tss =  Γ[1,1] / d
    tps = -Γ[3,1] / d
    tsp = -Γ[1,3] / d

    Rpp = abs2(rpp)
    Rss = abs2(rss)
    Rps = abs2(rps)
    Rsp = abs2(rsp)

    Tpp = abs2(tpp)
    Tss = abs2(tss)
    Tps = abs2(tps)
    Tsp = abs2(tsp)

    r = SVector(rpp, rps, rss, rsp)
    R = SVector(Rpp, Rss, Rsp, Rps)

    t = SVector(tpp, tps, tsp, tss)
    T = SVector(Tpp, Tss, Tsp, Tps)

    return r, R, t, T
end


"""
    calculate_tr(S::Poynting)

Calculate transmittance and reflectance from the Poynting vector struct,
which contains incident, transmitted, and reflected energy flux for both
p-polarized and s-polarized waves.

Returns `(Tpp, Tss, Rpp, Rss)`.

# Sign Convention
The reflected Poynting vector z-component is negative (pointing in -z direction),
so the negative sign in `Rpp = -S.refl_p[3] / S.in_p[3]` yields positive reflectance.
"""
function calculate_tr(S::Poynting)

    Tpp = S.out_p[3] / S.in_p[3]
    Tss = S.out_s[3] / S.in_s[3]

    # Reflected Poynting vector points in -z, so negate to get positive R
    Rpp = -S.refl_p[3] / S.in_p[3]
    Rss = -S.refl_s[3] / S.in_s[3]

    return Tpp, Tss, Rpp, Rss
end


"""
    transfer(λ, layers; θ=0.0, μ=1.0, validate=false, basis=:linear, method=:eig)

Calculate the transmittance and reflectance of a layered structure.

With the default `basis=:linear`, returns a [`TransferResult`](@ref) with fields
`Tpp`, `Tss`, `Rpp`, `Rss` (and cross-polarization terms `Tps`, `Tsp`, `Rps`, `Rsp`).
With `basis=:circular`, returns a [`CircularTransferResult`](@ref).

# Reflectance and transmittance calculation

**Reflectance** uses ``R = |r|^2`` from the transfer matrix coefficients
(Passler & Paarmann 2017, Eq. 17). This is exact for transparent incident
media. For absorbing incident media,
``|r|^2`` is not a true energy ratio — Poynting vectors become non-additive
due to interference cross-terms between incident and reflected waves
(Ortiz & Mochan 2005, JOSA A 22, 2827). Proper treatment of that case
requires the `power_entering` formalism (Byrnes 2016, arXiv:1603.02720),
which is not yet implemented here.

**Transmittance** uses Poynting vectors (energy flux ratio ``S_out / S_in``)
rather than ``|t|^2``, because the transmitted wave propagates in a different
medium than the incident wave. As noted in the 2019 erratum (JOSAB 36, 3246):
``T ≠ |t|^2`` in general; only when the substrate is vacuum does ``T = |t|^2``.

# Arguments
- `λ`: Wavelength in μm (must match units used for layer thicknesses)
- `layers`: Vector of `Layer` objects representing the stack
- `θ`: Angle of incidence in radians (default: 0.0, normal incidence)
- `μ`: Relative magnetic permeability (default: 1.0, non-magnetic)
- `validate`: Check energy conservation R + T ≈ 1 for non-absorbing media (default: false)
- `basis`: Output polarization basis — `:linear` (default) or `:circular`
- `method`: propagation backend — `:eig` (eigenmode) or `:exp` (matrix exponential, see [`layer_transfer_exp`](@ref))

# Polarization basis
- `basis=:linear` (default) returns a [`TransferResult`](@ref) in the linear p/s basis.
- `basis=:circular` returns a [`CircularTransferResult`](@ref) in the right/left
  circular basis. R/L are fixed-lab-frame helicities under this package's
  `exp(-iωt)` convention. The Jones matrices transform as `r_circ = C⁻¹ r_lin C`
  with `C = (1/√2)[1 1; i -i]` (columns ordered L, R); the helicity flip on
  reflection is encoded automatically in the opposite signs of `rpp`/`rss`.

  Circular **reflectance** is `|r_circ|²` (a true energy ratio under the same
  condition as linear `R`: transparent incident medium, cf. issue #72). Circular
  **transmittance** is `N·|t_circ|²` with a single Poynting normalization scalar
  `N = (Tpp+Tss)/(|tpp|²+|tss|²)`; this is exact and energy-conserving for an
  isotropic substrate (and reduces to `|t_circ|²` for a vacuum/index-matched
  substrate). For an **anisotropic substrate** the single scalar `N` is only
  approximate: the two transmitted eigenmodes carry different wavevectors, so
  their Poynting-to-`|t|²` conversion factors differ, and one polarization-
  independent scalar cannot capture both. (The linear basis avoids this by
  summing each transmitted eigenmode with its own wavevector in `poynting()`,
  so linear `T` — including the cross-pol `Tps/Tsp` — is per-mode exact.)
  Circular `T` for an anisotropic substrate is therefore not guaranteed to be
  energy-conserving to machine precision. `validate` applies to the linear basis only.

# Wave Propagation Convention
- Light propagates in the **+z direction** (from first layer toward last layer)
- The first and last layers are treated as semi-infinite media
- θ is measured from the surface normal (z-axis)

# Units
- Wavelength and thicknesses: μm by default. With `using Unitful`, `λ` may carry
  units — a length, or a wavenumber/frequency/photon-energy that maps to
  wavelength (e.g. `1550u"nm"`, `193u"THz"`, `0.8u"eV"`).
- Angle: radians
- Transmittance/Reflectance: dimensionless (0 to 1)

# Physics Validation
When `validate=true`, the function checks:
1. **Bounds**: 0 ≤ R, T ≤ 1 (catches NaN, negative values, numerical instability)
2. **Energy conservation**: R + T ≈ 1 for non-absorbing media (imag(n) < 1e-10)
3. **Absorption bound**: R + T ≤ 1 for absorbing media

Warnings are issued for any violations.
"""
function transfer(λ, layers; θ=0.0, μ=1.0, sheets=nothing, validate::Bool=false, basis::Symbol=:linear, method::Symbol=:eig)
    λ = _to_wavelength_um(λ)
    θ = _to_radians(θ)

    sd = sheets === nothing ? nothing : _sheets_dict(sheets)
    _validate_sheet_indices(sd, length(layers))
    Γ, S = _propagate(Val(method), λ, layers; θ=θ, μ=μ, sheets=sd)
    return _assemble(Val(basis), Γ, S, λ, layers, sd, validate)
end

function _assemble(::Val{:linear}, Γ, S, λ, layers, sd, validate)
    r, R, t, T = calculate_tr(Γ)
    Tpp, Tss, _, _ = calculate_tr(S)

    # Diagonal terms from matrix calculation; cross terms remapped from the
    # SVector packing (R = (Rpp,Rss,Rsp,Rps), T = (Tpp,Tss,Tsp,Tps)).
    Rpp = R[1]
    Rss = R[2]
    Rps = R[4]
    Rsp = R[3]
    Tps = T[4]
    Tsp = T[3]

    if validate
        _validate_physics(λ, layers, Tpp, Tss, Rpp, Rss; sheets=sd)
    end

    return TransferResult(Tpp, Tss, Tps, Tsp, Rpp, Rss, Rps, Rsp)
end

function _assemble(::Val{:circular}, Γ, S, λ, layers, sd, validate)
    r, _, t, _ = calculate_tr(Γ)
    Tpp, Tss, _, _ = calculate_tr(S)
    return _circular_result(r, t, Tpp, Tss)
end

_assemble(::Val{B}, Γ, S, λ, layers, sd, validate) where {B} =
    throw(ArgumentError("basis must be :linear or :circular, got :$B"))


"""
    _validate_physics(λ, layers, Tpp, Tss, Rpp, Rss; atol=1e-6, k_threshold=1e-10)

Validate physical constraints on R and T values:
1. Bounds check: 0 ≤ R, T ≤ 1 (always)
2. Energy conservation: R + T ≈ 1 (for non-absorbing media only)

Issues warnings if constraints are violated.

Internal function called by `transfer` when `validate=true`.
"""
function _validate_physics(λ, layers, Tpp, Tss, Rpp, Rss; sheets=nothing, atol=1e-6, k_threshold=1e-10)
    # Check for NaN values (indicates numerical failure)
    if any(isnan, (Tpp, Tss, Rpp, Rss))
        @warn "NaN detected in R/T values" Tpp Tss Rpp Rss
        return nothing
    end

    # Check physical bounds: 0 ≤ R, T ≤ 1
    if Tpp < 0 || Tpp > 1
        @warn "Transmittance Tpp out of bounds [0, 1]" Tpp
    end
    if Tss < 0 || Tss > 1
        @warn "Transmittance Tss out of bounds [0, 1]" Tss
    end
    if Rpp < 0 || Rpp > 1
        @warn "Reflectance Rpp out of bounds [0, 1]" Rpp
    end
    if Rss < 0 || Rss > 1
        @warn "Reflectance Rss out of bounds [0, 1]" Rss
    end

    # Check if all layers are non-absorbing
    layers_lossless = all(layers) do layer
        nx, ny, nz = get_refractive_indices(layer, λ)
        all(n -> abs(imag(n)) < k_threshold, (nx, ny, nz))
    end

    # Sheets are lossless when Re(σ) ≈ 0 (purely reactive). Any in-plane Re(σ) > 0 absorbs.
    sheets_lossless = sheets === nothing || all(values(sheets)) do sheet
        σ = sheet.conductivity(λ)
        all(c -> abs(real(c)) < k_threshold, (σ[1,1], σ[1,2], σ[2,1], σ[2,2]))
    end

    is_lossless = layers_lossless && sheets_lossless

    if is_lossless
        # Energy conservation: R + T = 1 for lossless media
        sum_p = Tpp + Rpp
        sum_s = Tss + Rss

        if !isapprox(sum_p, 1.0; atol=atol)
            @warn "Energy conservation violated for p-polarization" Tpp Rpp sum=sum_p expected=1.0 deviation=abs(sum_p - 1.0)
        end

        if !isapprox(sum_s, 1.0; atol=atol)
            @warn "Energy conservation violated for s-polarization" Tss Rss sum=sum_s expected=1.0 deviation=abs(sum_s - 1.0)
        end
    else
        # For lossy media: R + T ≤ 1
        sum_p = Tpp + Rpp
        sum_s = Tss + Rss

        if sum_p > 1.0 + atol
            @warn "Absorption violation: R + T > 1 for p-polarization" Tpp Rpp sum=sum_p
        end
        if sum_s > 1.0 + atol
            @warn "Absorption violation: R + T > 1 for s-polarization" Tss Rss sum=sum_s
        end
    end

    return nothing
end


function _sweep_spectra(outer_vals, inner_vals, ::Val{B}; threads::Bool=true, verbose::Bool=false, make_layers, angle_for, sheets=nothing) where {B}
    dims = (length(outer_vals), length(inner_vals))
    M1 = Array{Float64}(undef, dims)
    M2 = Array{Float64}(undef, dims)
    M3 = Array{Float64}(undef, dims)
    M4 = Array{Float64}(undef, dims)
    M5 = Array{Float64}(undef, dims)
    M6 = Array{Float64}(undef, dims)
    M7 = Array{Float64}(undef, dims)
    M8 = Array{Float64}(undef, dims)

    if verbose
        println("Threads: ", Threads.nthreads())
    end

    function compute_row(i)
        layers_i = make_layers(i)
        θ = angle_for(i)
        for j in eachindex(inner_vals)
            result = transfer(inner_vals[j], layers_i; θ=θ, sheets=sheets, basis=B)
            if B === :linear
                M1[i, j] = result.Tpp
                M2[i, j] = result.Tss
                M3[i, j] = result.Tps
                M4[i, j] = result.Tsp
                M5[i, j] = result.Rpp
                M6[i, j] = result.Rss
                M7[i, j] = result.Rps
                M8[i, j] = result.Rsp
            else
                M1[i, j] = result.Trr
                M2[i, j] = result.Tll
                M3[i, j] = result.Trl
                M4[i, j] = result.Tlr
                M5[i, j] = result.Rrr
                M6[i, j] = result.Rll
                M7[i, j] = result.Rrl
                M8[i, j] = result.Rlr
            end
        end
    end

    if threads
        Threads.@threads for i in eachindex(outer_vals)
            compute_row(i)
        end
    else
        for i in eachindex(outer_vals)
            compute_row(i)
        end
    end

    return B === :linear ?
        TransferResult(M1, M2, M3, M4, M5, M6, M7, M8) :
        CircularTransferResult(M1, M2, M3, M4, M5, M6, M7, M8)
end

"""
    sweep_angle(λs, θs, layers; threads=true, verbose=false, basis=:linear)

Calculate transmittance/reflectance spectra over wavelength and angle of incidence.

Returns a `TransferResult` with fields `Tpp`, `Tss`, `Rpp`, `Rss`, each a matrix
of size `(length(θs), length(λs))`.

# Arguments
- `λs`: Vector of wavelengths in μm
- `θs`: Vector of angles of incidence in radians
- `layers`: `AbstractVector{<:Layer}` representing the stack
- `threads`: Enable multithreading (default: true)
- `verbose`: Print thread count info (default: false)
- `basis`: `:linear` (default) → `TransferResult`; `:circular` → `CircularTransferResult` (see [`transfer`](@ref))

# Units
- Wavelengths: μm (micrometers) recommended
- Angles: radians
"""
function sweep_angle(λs, θs, layers; sheets=nothing, threads::Bool=true, verbose::Bool=false, basis::Symbol=:linear)
    λs = _to_wavelength_um.(λs)
    θs = _to_radians.(θs)
    sd = sheets === nothing ? nothing : _sheets_dict(sheets)
    _validate_sheet_indices(sd, length(layers))
    return _sweep_spectra(θs, λs, Val(basis); threads=threads, verbose=verbose,
        make_layers = _ -> layers,
        angle_for = i -> θs[i],
        sheets = sd)
end


"""
    sweep_thickness(λs, ts, layers, t_index; θ=0.0, threads=true, verbose=false, basis=:linear)

Sweep the thickness of a specific layer and calculate transmittance/reflectance spectra.

Returns a `TransferResult` with fields `Tpp`, `Tss`, `Rpp`, `Rss`, each a matrix
of size `(length(ts), length(λs))`.

# Arguments
- `λs`: Vector of wavelengths in μm
- `ts`: Vector of thicknesses in μm to sweep
- `layers`: `AbstractVector{<:Layer}` representing the stack
- `t_index`: Index of the layer (1-based) whose thickness to vary
- `θ`: Angle of incidence in radians (default: 0.0, normal incidence)
- `threads`: Enable multithreading (default: true)
- `verbose`: Print thread count info (default: false)
- `basis`: `:linear` (default) → `TransferResult`; `:circular` → `CircularTransferResult` (see [`transfer`](@ref))

# Units
- Wavelengths and thicknesses: μm (micrometers) recommended
- Angle: radians
"""
function sweep_thickness(λs, ts, layers, t_index::Int; θ=0.0, sheets=nothing, threads::Bool=true, verbose::Bool=false, basis::Symbol=:linear)
    λs = _to_wavelength_um.(λs)
    ts = _to_um.(ts)
    θ = _to_radians(θ)
    sd = sheets === nothing ? nothing : _sheets_dict(sheets)
    _validate_sheet_indices(sd, length(layers))
    dispersion_func = layers[t_index].dispersion
    layers_base = collect(layers)

    mu_orig = layers[t_index].mu
    return _sweep_spectra(ts, λs, Val(basis); threads=threads, verbose=verbose,
        make_layers = i -> begin
            layers_i = copy(layers_base)
            layers_i[t_index] = Layer(dispersion_func, mu_orig, ts[i])
            layers_i
        end,
        angle_for = _ -> θ,
        sheets = sd)
end

@deprecate angle_resolved(λs, θs, layers; kwargs...) sweep_angle(λs, θs, layers; kwargs...)
@deprecate tune_thickness(λs, ts, layers, t_index::Int, θ=0.0; kwargs...) sweep_thickness(λs, ts, layers, t_index; θ=θ, kwargs...)
@deprecate calculate_tr(λ, layers; kwargs...) transfer(λ, layers; kwargs...)
@deprecate electric_field(λ, layers; kwargs...) efield(λ, layers; kwargs...)


# Shared core for efield/hfield: runs _propagate_full once, performs the backward
# mode-coefficient recursion (with sheet injection), samples the z-grid, and returns
# everything both wrappers need. E and H differ only in the final per-z reconstruction.
function _field(λ, layers; θ=0.0, μ=1.0, dz=0.001, sheets=nothing)
    λ = _to_wavelength_um(λ)
    dz = _to_um(dz)
    θ = _to_radians(θ)

    sd = sheets === nothing ? nothing : _sheets_dict(sheets)
    _validate_sheet_indices(sd, length(layers))
    no_sheets = sd === nothing || isempty(sd)

    Γ, S, Ds, Ps, γs, qs = _propagate_full(λ, layers; θ=θ, μ=μ, sheets=sd)
    r, R, t, T = calculate_tr(Γ)

    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ξ = √(dielectric_constant(nx_in)) * sin(θ)

    first_layer = layers[1]
    last_layer = layers[end]
    nlay = length(layers)

    Eplus_p = zeros(ComplexF64, nlay, 4)
    Eminus_p = zeros(ComplexF64, nlay, 4)
    Eplus_s = zeros(ComplexF64, nlay, 4)
    Eminus_s = zeros(ComplexF64, nlay, 4)

    Eplus_p[end, :] = [t[1], t[2], 0, 0]
    Eplus_s[end, :] = [t[3], t[4], 0, 0]

    P_f = Ps[end]
    Eminus_p[end, :] = P_f(last_layer.thickness) \ Eplus_p[end, :]
    Eminus_s[end, :] = P_f(last_layer.thickness) \ Eplus_s[end, :]

    D_i = Ds[end]
    for l in reverse(eachindex(layers))
        if l >= 2
            layer = layers[l - 1]
            P_prev = Ps[l - 1]
            if no_sheets || !haskey(sd, l - 1)
                L_i = Ds[l - 1] \ D_i
            else
                L_i = Ds[l - 1] \ (sheet_matrix(sd[l - 1], λ) * D_i)
            end
            Eminus_p[l - 1, :] = L_i * Eplus_p[l, :]
            Eminus_s[l - 1, :] = L_i * Eplus_s[l, :]
            Eplus_p[l - 1, :] = P_prev(layer.thickness) * Eminus_p[l - 1, :]
            Eplus_s[l - 1, :] = P_prev(layer.thickness) * Eminus_s[l - 1, :]
            D_i = Ds[l - 1]
        end
    end

    interface_positions, total_thickness = find_bounds(layers)
    interface_positions .-= first_layer.thickness
    zs = range(-first_layer.thickness, interface_positions[end], step=dz)

    nz = length(zs)
    amp_p = zeros(ComplexF64, 4, nz)
    amp_s = zeros(ComplexF64, 4, nz)
    layer_of_z = Vector{Int}(undef, nz)

    i = 1
    for (j, z) in enumerate(zs)
        if i < nlay && z > interface_positions[i]
            i += 1
        end
        P_i = Ps[i]
        amp_p[:, j] = P_i(-(z - interface_positions[i])) * Eminus_p[i, :]
        amp_s[:, j] = P_i(-(z - interface_positions[i])) * Eminus_s[i, :]
        layer_of_z[j] = i
    end

    μ_mat = SMatrix{3,3,ComplexF64}(μ * I)
    μs = [ismagnetic(L) ? get_permeability(L, λ) : μ_mat for L in layers]

    return (; zs, boundaries = interface_positions[1:end - 1],
              amp_p, amp_s, layer_of_z, γs, qs, ξ, μ, μs)
end


"""
    efield(λ, layers; θ=0.0, μ=1.0, dz=0.001, sheets=nothing)

Calculate the electric field profile throughout the layered structure.

Returns an `ElectricField` struct containing:
- `z`: Position coordinates along the structure
- `p`: Electric field components (Ex, Ey, Ez) for p-polarized incidence
- `s`: Electric field components (Ex, Ey, Ez) for s-polarized incidence
- `boundaries`: z-positions of layer interfaces

# Arguments
- `λ`: Wavelength in μm (must match units used for layer thicknesses)
- `layers`: Vector of `Layer` objects representing the stack
- `θ`: Angle of incidence in radians (default: 0.0, normal incidence)
- `μ`: Relative magnetic permeability (default: 1.0, non-magnetic)
- `dz`: Spatial step size in μm for field sampling (default: 0.001)
- `sheets`: Optional conductive sheets at interfaces (see [`Sheet`](@ref) and
  [`transfer`](@ref)); keyed by the index of the layer above each interface.

# Wave Propagation Convention
- Light propagates in the **+z direction** (from first layer toward last layer)
- z = 0 is at the first interface (between layer 1 and layer 2)
- Negative z values are inside the first (incident) layer
- θ is measured from the surface normal (z-axis)

# Units
- All lengths (λ, thickness, dz, z): μm by default; with `using Unitful`, `λ`
  and `dz` may carry units (e.g. `efield(1.55u"μm", layers; dz=1u"nm")`).
- Angle: radians
- Electric field: arbitrary units (normalized to incident field)
"""
function efield(λ, layers; θ=0.0, μ=1.0, dz=0.001, sheets=nothing)

    F = _field(λ, layers; θ=θ, μ=μ, dz=dz, sheets=sheets)
    nz = length(F.zs)
    p = zeros(ComplexF64, 3, nz)
    s = zeros(ComplexF64, 3, nz)

    for j in 1:nz
        γ = F.γs[F.layer_of_z[j]]
        ap = view(F.amp_p, :, j)
        as = view(F.amp_s, :, j)
        @views p[:, j] = ap[1] * γ[1, :] + ap[2] * γ[2, :] + ap[3] * γ[3, :] + ap[4] * γ[4, :]
        @views s[:, j] = as[1] * γ[1, :] + as[2] * γ[2, :] + as[3] * γ[3, :] + as[4] * γ[4, :]
    end

    return ElectricField(F.zs, p, s, F.boundaries)
end


# H eigenvectors per mode from the E eigenvectors γ and eigenvalues q:
# H_m = μ⁻¹(k̄×E)_m = μ⁻¹(-q γ₂, q γ₁ - ξ γ₃, ξ γ₂) = (Hx, Hy, Hz).
# Applies the full μ⁻¹ tensor. Rows 2,1 match dynamical_matrix rows 3,4
# (H_y and -Hx); row 3 (Hz) is (k×E)_z = ξ E_y.
function _h_eigvecs(γ, q, ξ, μ::AbstractMatrix)
    μinv = inv(SMatrix{3,3,ComplexF64}(μ))
    η = @MMatrix zeros(ComplexF64, 4, 3)
    for m in 1:4
        E = SVector{3,ComplexF64}(γ[m, 1], γ[m, 2], γ[m, 3])
        H = μinv * (_kcross(ξ, q[m]) * E)
        η[m, 1] = H[1]; η[m, 2] = H[2]; η[m, 3] = H[3]
    end
    return SMatrix(η)
end

"""
    hfield(λ, layers; θ=0.0, μ=1.0, dz=0.001, sheets=nothing)

Calculate the magnetic-field profile through the structure, returning a
[`MagneticField`](@ref). Shares the sampling grid with [`efield`](@ref) (same `dz`
and arguments), so E and H can be overlaid; the in-plane H discontinuity at a
conductive sheet equals the surface current `ẑ × (σ_s E∥)`.

# Units / normalization
H is returned in impedance-normalized units `H̃ = Z₀ H_SI` (`Z₀ = √(μ₀/ε₀)`), so
`|E| ~ |H̃|` for a plane wave. Arguments and conventions match [`efield`](@ref).
"""
function hfield(λ, layers; θ=0.0, μ=1.0, dz=0.001, sheets=nothing)

    F = _field(λ, layers; θ=θ, μ=μ, dz=dz, sheets=sheets)
    nz = length(F.zs)
    p = zeros(ComplexF64, 3, nz)
    s = zeros(ComplexF64, 3, nz)

    # η depends only on the layer index (via γ and q), not on the z-sample, so
    # build it once per layer and index by layer rather than rebuilding per z.
    ηs = [_h_eigvecs(F.γs[li], F.qs[li], F.ξ, F.μs[li]) for li in eachindex(F.γs)]

    for j in 1:nz
        η = ηs[F.layer_of_z[j]]
        ap = view(F.amp_p, :, j)
        as = view(F.amp_s, :, j)
        @views p[:, j] = ap[1] * η[1, :] + ap[2] * η[2, :] + ap[3] * η[3, :] + ap[4] * η[4, :]
        @views s[:, j] = as[1] * η[1, :] + as[2] * η[2, :] + as[3] * η[3, :] + as[4] * η[4, :]
    end

    return MagneticField(F.zs, p, s, F.boundaries)
end
