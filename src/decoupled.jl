# 2×2 fast path for stacks in which s and p light never mix.
#
# Berreman's 4×4 Δ splits into a p block acting on (Eₓ, Hᵧ) and an s block acting
# on (Eᵧ, −Hₓ) whenever y (the axis normal to the plane of incidence) is a
# principal axis of every ε and μ in the stack and every sheet conductivity is
# diagonal in (x, y). A product of block-diagonal matrices stays block-diagonal,
# so the whole system matrix is diag(Γp, Γs) and the cross-polarized coefficients
# vanish. Each 2×2 block has a closed-form propagator, which removes the per-layer
# matrix exponential and the ambient/substrate eigen-solves of the 4×4 path.
#
# This file covers the cases decidable from the layer structure alone:
#   • ambient and substrate isotropic and non-magnetic;
#   • interior layers isotropic, or anisotropic with zero Euler angles
#     (principal axes on the lab axes), with a diagonal μ tensor if magnetic;
#   • sheets with σxy = σyx = 0.
# Rotated layers, gyrotropic μ and Hall-type sheets go to the 4×4 path.
#
# The ambient and substrate use the package's own `calculate_E_modes` and
# `dynamical_matrix` with closed-form mode wavevectors, so `M_sys` and the
# Poynting struct are expressed in exactly the basis the 4×4 path uses: every
# consumer of `(M_sys, S)` (r/t, R/T, circular basis, amplitude coefficients,
# the system matrix itself) gets the same numbers to rounding.

# Relative size below which an off-block tensor element counts as zero. Rotation
# by Euler angles such as π/2 leaves ~1e-16 residue; a genuine coupling of
# relative size δ moves co-polarized outputs only at O(δ²).
const _DECOUPLE_RTOL = 1e-12

_isotropic_boundary(layer) = !isanisotropic(layer) && !ismagnetic(layer)

"""
    _decoupled_structure(layers)

Structural pre-check for the 2×2 path: isotropic, non-magnetic ambient and
substrate, and no rotated interior layer. Per-point conditions (diagonal μ,
diagonal sheet conductivity) are checked inside [`_propagate_core_2x2`](@ref).
"""
function _decoupled_structure(layers)
    N = length(layers)
    N ≥ 2 || return false
    (_isotropic_boundary(layers[1]) && _isotropic_boundary(layers[N])) || return false
    for i in 2:N-1
        isrotated(layers[i]) && return false
    end
    return true
end

# Diagonal permeability (μx, μy, μz) of an interior layer at λ, or `nothing` if the
# layer's μ tensor has off-diagonal elements (gyrotropic media couple s and p).
@inline function _diagonal_mu(layer, λ, μ_i)
    if !ismagnetic(layer)
        m = ComplexF64(μ_i)
        return (m, m, m)
    end
    M = get_permeability(layer, λ)
    off = max(abs(M[1,2]), abs(M[2,1]), abs(M[1,3]), abs(M[3,1]), abs(M[2,3]), abs(M[3,2]))
    off <= _DECOUPLE_RTOL * maximum(abs, M) || return nothing
    return (M[1,1], M[2,2], M[3,3])
end

# exp(-i x A) for A = [0 a; b 0]: A² = ab·I, so with q = √(ab)
#   exp(-i x A) = cos(xq) I - i sin(xq)/q A,
# which is even in q (no branch choice) and finite as q → 0.
@inline function _cos_sinc(x, q)
    c = cos(x * q)
    s = iszero(q) ? complex(0.0, -x) : -im * sin(x * q) / q
    return c, s
end
@inline _block_from(a, b, c, s) = SMatrix{2,2,ComplexF64}(c, s * b, s * a, c)   # [c s·a; s·b c]
@inline function _block_propagator(a, b, x)
    c, s = _cos_sinc(x, sqrt(a * b))
    return _block_from(a, b, c, s)
end

# Forward branch of a mode wavevector, following `calculate_q`'s sorting rule:
# effectively-real modes go forward with Re q ≥ 0, others with Im q ≥ 0.
@inline function _forward_branch(q)
    if abs(imag(q)) < 1e-10 * max(abs(q), 1.0)
        return real(q) >= 0 ? q : -q
    end
    return imag(q) >= 0 ? q : -q
end

# Closed-form replacement for `layer_matrices` on an isotropic, non-magnetic
# semi-infinite medium: the mode wavevectors are ±√(με − k∥²), twice degenerate,
# and the E modes and dynamical matrix come from the package's scalar-μ builders.
function _isotropic_boundary_modes(layer, λ, k_par, μ_i)
    n, _, _ = get_refractive_indices(layer, λ)
    ε = dielectric_constant(n)
    qz = _forward_branch(sqrt(complex(μ_i * ε - k_par^2)))
    q = SVector{4,ComplexF64}(qz, qz, -qz, -qz)
    E_modes = calculate_E_modes(k_par, q, dielectric_tensor(ε, ε, ε), μ_i)
    D = dynamical_matrix(k_par, q, E_modes, μ_i)
    return D, E_modes, q
end

# p block (D-basis rows/cols 1, 3) and s block (rows/cols 2, 4) of a sheet matrix,
# or `nothing` if the sheet couples them (σxy or σyx nonzero).
@inline function _sheet_blocks(sheet, λ)
    G = sheet_matrix(sheet, λ)
    off = max(abs(G[3,2]), abs(G[4,1]))
    off <= _DECOUPLE_RTOL * max(abs(G[3,1]), abs(G[4,2])) || return nothing
    one_ = one(ComplexF64); z = zero(ComplexF64)
    return SMatrix{2,2,ComplexF64}(one_, G[3,1], z, one_), SMatrix{2,2,ComplexF64}(one_, G[4,2], z, one_)
end

# p block (rows/cols 1, 3) and s block (rows/cols 2, 4) of a 4×4 matrix in the
# dynamical-matrix basis (Eₓ, Eᵧ, Hᵧ, −Hₓ).
@inline _pblock(D) = SMatrix{2,2,ComplexF64}(D[1,1], D[3,1], D[1,3], D[3,3])
@inline _sblock(D) = SMatrix{2,2,ComplexF64}(D[2,2], D[4,2], D[2,4], D[4,4])

const _Block2 = SMatrix{2,2,ComplexF64,4}
const _BlockResult = Tuple{Bool,_Block2,_Block2}

# p and s propagators of one interior layer, or ok = false if its μ couples s and p.
# Berreman's blocks for diagonal ε = (εx, εy, εz), μ = (μx, μy, μz):
#   p: [0  μy − k∥²/εz; εx  0]      s: [0  μx; εy − k∥²/μz  0]
function _interior_blocks(layer, λ, k2, k0, μ_i)::_BlockResult
    mu = _diagonal_mu(layer, λ, μ_i)
    mu === nothing && return (false, zero(_Block2), zero(_Block2))
    μx, μy, μz = mu
    nx, ny, nz = get_refractive_indices(layer, λ)
    εx = dielectric_constant(nx); εy = dielectric_constant(ny); εz = dielectric_constant(nz)
    x = k0 * layer.thickness
    ap, bp = μy - k2 / εz, εx
    as, bs = μx, εy - k2 / μz
    if !isanisotropic(layer) && μx == μy == μz
        # Isotropic: both blocks share q² = με − k∥², so one cos/sin serves both.
        c, s = _cos_sinc(x, sqrt(ap * bp))
        return (true, _block_from(ap, bp, c, s), _block_from(as, bs, c, s))
    end
    return (true, _block_propagator(ap, bp, x), _block_propagator(as, bs, x))
end

"""
    _propagate_core_2x2(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

2×2 propagation core for decoupled stacks. Returns `(M_sys, S)` in the same basis
as [`_propagate_core_exp`](@ref), or `nothing` when a per-point check fails (a
magnetic interior layer with off-diagonal μ, or a sheet with σxy ≠ 0), in which
case the caller falls back to the 4×4 path. Assumes `_decoupled_structure(layers)`.

Per interior layer with diagonal ε = (εx, εy, εz) and μ = (μx, μy, μz), Berreman's
blocks are p: [0  μy − k∥²/εz; εx  0] and s: [0  μx; εy − k∥²/μz  0].
"""
function _propagate_core_2x2(λ, layers; θ=0.0, μ=1.0, sheets=nothing)
    N = length(layers)
    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ε_0in = dielectric_constant(nx_in)
    k_par = √(ε_0in) * sin(θ)
    k2 = k_par^2
    k0 = (2π * c_0 / λ) / c_0
    no_sheets = sheets === nothing || isempty(sheets)

    Bp = SMatrix{2,2,ComplexF64}(I)
    Bs = SMatrix{2,2,ComplexF64}(I)
    if !no_sheets && haskey(sheets, 1)
        blocks = _sheet_blocks(sheets[1], λ)
        blocks === nothing && return nothing
        Bp = Bp * blocks[1]; Bs = Bs * blocks[2]
    end
    for i in 2:N-1
        # Function barrier: a stack mixing dispersion types is a Vector of abstract
        # Layer, so dispatch once per layer rather than once per field access.
        ok, P, Sb = _interior_blocks(layers[i], λ, k2, k0, μ)::_BlockResult
        ok || return nothing
        Bp = Bp * P
        Bs = Bs * Sb
        if !no_sheets && haskey(sheets, i)
            blocks = _sheet_blocks(sheets[i], λ)
            blocks === nothing && return nothing
            Bp = Bp * blocks[1]; Bs = Bs * blocks[2]
        end
    end

    D_1, E_modes_first, q_first = _isotropic_boundary_modes(layers[1], λ, k_par, μ)
    D_N, E_modes_last,  q_last  = _isotropic_boundary_modes(layers[N], λ, k_par, μ)

    # The isotropic dynamical matrices are block-diagonal too (p modes 1, 3 carry no
    # Eᵧ or Hₓ; s modes 2, 4 no Eₓ, E_z or Hᵧ), so D₁⁻¹·B·D_N splits into two 2×2
    # products. In the swapped basis of `_propagate_core_exp` the p block lands in
    # rows/columns 1:2 and the s block in 3:4.
    Mp = inv(_pblock(D_1)) * Bp * _pblock(D_N)
    Ms = inv(_sblock(D_1)) * Bs * _sblock(D_N)
    z = zero(ComplexF64)
    M_sys = @SMatrix [Mp[1,1] Mp[1,2] z z;
                      Mp[2,1] Mp[2,2] z z;
                      z z Ms[1,1] Ms[1,2];
                      z z Ms[2,1] Ms[2,2]]
    r, R, t, T = calculate_tr(M_sys)
    μ_mat = SMatrix{3,3,ComplexF64}(μ * I)
    S = poynting(k_par, q_first, q_last, E_modes_first, E_modes_last, t, r, μ_mat, μ_mat)
    return M_sys, S
end

"""
    isdecoupled(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Return `true` if `transfer(...; method=:auto)` would use the 2×2 decoupled path for
this point, i.e. s and p polarizations propagate independently through the stack
and the cross-polarized coefficients are exactly zero.
"""
function isdecoupled(λ, layers; θ=0.0, μ=1.0, sheets=nothing)
    _decoupled_structure(layers) || return false
    sd = sheets === nothing ? nothing : _sheets_dict(sheets)
    return _propagate_core_2x2(_to_wavelength_um(λ), layers; θ=_to_radians(θ), μ=μ, sheets=sd) !== nothing
end
