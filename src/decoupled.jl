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
# The ambient and substrate must be isotropic and non-magnetic. Interior layers are
# handled two ways:
#   • isotropic, or anisotropic with zero Euler angles and a diagonal μ: closed-form
#     blocks built straight from (εx, εy, εz) and (μx, μy, μz);
#   • anything else (rotated ε, full μ tensor): the layer's Berreman Δ is built and
#     its off-diagonal blocks are checked numerically at this point. If they vanish,
#     Δ's two diagonal blocks are exponentiated in closed form. This covers an optic
#     axis tilted in the plane of incidence (εxz ≠ 0), biaxial crystals with a
#     principal axis on y, μ with gyration axis y, and, at normal incidence only,
#     layers whose coupling enters through ξ alone (εyz, μyz).
# Sheets must have σxy = σyx = 0. Everything else goes to the 4×4 path.
#
# The ambient and substrate use the package's own `calculate_E_modes` and
# `dynamical_matrix` with closed-form mode wavevectors, so `M_sys` and the
# Poynting struct are expressed in exactly the basis the 4×4 path uses: every
# consumer of `(M_sys, S)` (r/t, R/T, circular basis, amplitude coefficients,
# the system matrix itself) gets the same numbers to rounding.

# Relative size below which an off-block tensor (or Δ) element counts as zero. Rotation
# by Euler angles such as π/2 leaves ~1e-16 residue; a genuine coupling of
# relative size δ moves co-polarized outputs only at O(δ²).
const _DECOUPLE_RTOL = 1e-12

_isotropic_boundary(layer) = !isanisotropic(layer) && !ismagnetic(layer)

"""
    _decoupled_structure(layers)

Structural pre-check for the 2×2 path: an isotropic, non-magnetic ambient and
substrate. Interior layers and sheets are checked per point inside
[`_propagate_core_2x2`](@ref), since whether a rotated or magnetic layer couples s
and p can depend on λ and θ.
"""
function _decoupled_structure(layers)
    N = length(layers)
    N ≥ 2 || return false
    return _isotropic_boundary(layers[1]) && _isotropic_boundary(layers[N])
end

# Diagonal permeability (μx, μy, μz) of an interior layer at λ, or `nothing` if the
# layer's μ tensor has off-diagonal elements (gyrotropic media couple s and p).
@inline function _diagonal_mu(layer, λ, μ_i)
    if !ismagnetic(layer)
        m = complex(float(μ_i))
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
    s = iszero(q) ? complex(zero(x), -x) : -im * sin(x * q) / q
    return c, s
end
# The element type C of every 2×2 block is fixed by the caller (`_stack_ctype`), so
# the per-layer function barrier below has a concrete return type even when some
# layers carry dual numbers and others do not.
@inline _block_from(::Type{C}, a, b, c, s) where {C} = SMatrix{2,2,C}(c, s * b, s * a, c)   # [c s·a; s·b c]
@inline function _block_propagator(::Type{C}, a, b, x) where {C}
    c, s = _cos_sinc(x, sqrt(a * b))
    return _block_from(C, a, b, c, s)
end

# exp(-i x A) for a general 2×2 A. With t = tr(A)/2 and B = A − tI (traceless),
# B² = q² I where q² = B₁₁² + B₁₂B₂₁, so exp(-i x A) = e^{-ixt} (cos(xq) I − i sin(xq)/q B).
@inline function _expm2(::Type{C}, A::SMatrix{2,2}, x) where {C}
    t = (A[1,1] + A[2,2]) / 2
    b11 = A[1,1] - t
    c, s = _cos_sinc(x, sqrt(b11^2 + A[1,2] * A[2,1]))
    e = exp(-im * x * t)
    return SMatrix{2,2,C}(e * (c + s * b11), e * s * A[2,1], e * s * A[1,2], e * (c - s * b11))
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
    return _isotropic_modes(dielectric_constant(n), k_par, μ_i)
end
#
# This is `calculate_E_modes` and the scalar-μ `dynamical_matrix` specialized to
# ε = ε·I with q₁ = q₂ = q_z, q₃ = q₄ = −q_z: the p modes are (1, 0, ∓k∥q_z/(με − k∥²))
# up to normalization (with the backward p mode's Eₓ = −1) and the s modes are ŷ.
function _isotropic_modes(ε, k_par, μ_i)
    qz = _forward_branch(sqrt(complex(μ_i * ε - k_par^2)))
    C = complex(float(promote_type(typeof(qz), typeof(ε), typeof(k_par), typeof(μ_i))))
    q = SVector{4,C}(qz, qz, -qz, -qz)
    o = one(C); z = zero(C)
    denom = μ_i * ε - k_par^2
    singular = abs(denom) < eps(Float64)
    e13 = singular ? z : (-k_par * qz) / denom
    e33 = singular ? z : (k_par * (-qz)) / denom
    v1 = SVector{3,C}(o, z, e13); v1 = v1 / norm(v1)
    v3 = SVector{3,C}(-o, z, e33); v3 = v3 / norm(v3)
    E_modes = SMatrix{4,3,C}(v1[1], z, v3[1], z,
                             z,     o, z,     o,
                             v1[3], z, v3[3], z)
    A3 = (E_modes[:, 1] .* q .- k_par * E_modes[:, 3]) ./ μ_i
    A4 = E_modes[:, 2] .* q ./ μ_i
    D = SMatrix{4,4,C}(E_modes[1,1], E_modes[1,2], A3[1], A4[1],
                       E_modes[2,1], E_modes[2,2], A3[2], A4[2],
                       E_modes[3,1], E_modes[3,2], A3[3], A4[3],
                       E_modes[4,1], E_modes[4,2], A3[4], A4[4])
    return D, E_modes, q
end

# What the 2×2 path needs of an isotropic, non-magnetic semi-infinite medium: the p
# and s blocks of its dynamical matrix and the Poynting vectors of its four unit
# modes (forward p, forward s, backward p, backward s), all in the basis of
# `_isotropic_modes`, so `(M_sys, S)` match the 4×4 path to rounding.
struct _IsoMedium{C<:Complex,R<:Real}
    Dp::SMatrix{2,2,C,4}
    Ds::SMatrix{2,2,C,4}
    Sp_f::SVector{3,R}
    Ss_f::SVector{3,R}
    Sp_b::SVector{3,R}
    Ss_b::SVector{3,R}
end

# Time-averaged Poynting vectors ½ Re(E × H*), H = (k × E)/μ with k = (k∥, 0, q)/c₀
# (the expression used by `poynting`), written out for the two mode shapes:
# a p mode E = (a, 0, b) has H = (0, qa − k∥b, 0)/(c₀μ); the s mode E = ŷ has
# H = (−q, 0, k∥)/(c₀μ).
@inline function _unit_poynting_p(a, b, k_par, q, μ)
    h = conj((q * a - k_par * b) / (c_0 * μ))
    return SVector(-real(b * h), zero(real(h)), real(a * h)) / 2
end
@inline function _unit_poynting_s(k_par, q, μ)
    return SVector(real(k_par / μ), zero(real(q)), real(q / μ)) / (2c_0)
end

function _iso_medium(::Type{C}, ε, k_par, μ) where {C}
    qz = _forward_branch(sqrt(complex(μ * ε - k_par^2)))
    o = one(C); z = zero(C)
    denom = μ * ε - k_par^2
    singular = abs(denom) < eps(Float64)
    e13 = singular ? z : (-k_par * qz) / denom
    e33 = singular ? z : (k_par * (-qz)) / denom
    # unit-normalized p modes (a, 0, b): forward (1, 0, e13), backward (−1, 0, e33)
    n1 = sqrt(1 + abs2(e13)); a1 = o / n1; b1 = e13 / n1
    n3 = sqrt(1 + abs2(e33)); a3 = -o / n3; b3 = e33 / n3
    Dp = SMatrix{2,2,C}(a1, (a1 * qz - k_par * b1) / μ,
                        a3, (a3 * (-qz) - k_par * b3) / μ)
    Ds = SMatrix{2,2,C}(o, qz / μ, o, -qz / μ)
    return _IsoMedium(Dp, Ds, _unit_poynting_p(a1, b1, k_par, qz, μ), _unit_poynting_s(k_par, qz, μ),
                      _unit_poynting_p(a3, b3, k_par, -qz, μ), _unit_poynting_s(k_par, -qz, μ))
end

# Poynting struct of a decoupled stack from its boundary media and r/t coefficients
# (r = (rpp, rps, rss, rsp), t = (tpp, tps, tsp, tss) as from `calculate_tr`). Mode
# fluxes scale with |amplitude|², and the cross-polarized amplitudes are exactly
# zero here, so every entry is |coefficient|² times a unit-mode vector, as `poynting`
# would compute it.
function _poynting_decoupled(a::_IsoMedium, b::_IsoMedium, t, r)
    out_pp = abs2(t[1]) * b.Sp_f
    out_ps = abs2(t[2]) * b.Ss_f
    out_sp = abs2(t[3]) * b.Sp_f
    out_ss = abs2(t[4]) * b.Ss_f
    refl_p = abs2(r[1]) * a.Sp_b + abs2(r[2]) * a.Ss_b
    refl_s = abs2(r[3]) * a.Sp_b + abs2(r[4]) * a.Ss_b
    return Poynting(out_pp + out_ps, a.Sp_f, out_sp + out_ss, a.Ss_f, refl_p, refl_s,
                    out_pp, out_ps, out_sp, out_ss, b.Sp_f[3], b.Ss_f[3])
end

# p block (D-basis rows/cols 1, 3) and s block (rows/cols 2, 4) of a sheet matrix,
# or `nothing` if the sheet couples them (σxy or σyx nonzero).
@inline function _sheet_blocks(::Type{C}, sheet, λ) where {C}
    G = sheet_matrix(sheet, λ)
    off = max(abs(G[3,2]), abs(G[4,1]))
    off <= _DECOUPLE_RTOL * max(abs(G[3,1]), abs(G[4,2])) || return nothing
    one_ = one(C); z = zero(C)
    return SMatrix{2,2,C}(one_, G[3,1], z, one_), SMatrix{2,2,C}(one_, G[4,2], z, one_)
end

# p block (rows/cols 1, 3) and s block (rows/cols 2, 4) of a 4×4 matrix in the
# dynamical-matrix basis (Eₓ, Eᵧ, Hᵧ, −Hₓ).
@inline _pblock(D) = SMatrix{2,2}(D[1,1], D[3,1], D[1,3], D[3,3])
@inline _sblock(D) = SMatrix{2,2}(D[2,2], D[4,2], D[2,4], D[4,4])

# Status of a per-layer 2×2 step: decoupled, couples s and p (fall back to 4×4), or
# the element type C is too narrow for the layer's indices (recompute, see `_fits`).
const _DECOUPLED = 0x00
const _COUPLED = 0x01
const _NARROW = 0x02
const _BlockResult{C} = Tuple{UInt8,SMatrix{2,2,C,4},SMatrix{2,2,C,4}}

# Berreman Δ (order Eₓ, Hᵧ, Eᵧ, −Hₓ) of a rotated or full-μ layer, and whether its
# off-diagonal blocks vanish at this point.
function _tensor_delta(layer, λ, k_par, μ_i, ns=get_refractive_indices(layer, λ))
    ε = _layer_epsilon(layer, ns...)
    μ = ismagnetic(layer) ? get_permeability(layer, λ) : permeability_tensor(μ_i, μ_i, μ_i)
    C = construct_constitutive(ε, μ)
    Δ = SMatrix{4,4}(construct_Δ(k_par, C, construct_a(k_par, C)))
    off = max(abs(Δ[1,3]), abs(Δ[1,4]), abs(Δ[2,3]), abs(Δ[2,4]),
              abs(Δ[3,1]), abs(Δ[3,2]), abs(Δ[4,1]), abs(Δ[4,2]))
    return off <= _DECOUPLE_RTOL * maximum(abs, Δ), Δ
end

# p and s propagators of a rotated or full-μ layer from the diagonal blocks of its
# Berreman Δ, or ok = false if the off-diagonal blocks do not vanish here.
function _tensor_blocks(::Type{C}, layer, λ, k_par, k0, μ_i, ns)::_BlockResult{C} where {C}
    ok, Δ = _tensor_delta(layer, λ, k_par, μ_i, ns)
    Z = zero(SMatrix{2,2,C,4})
    ok || return (_COUPLED, Z, Z)
    x = k0 * layer.thickness
    Δp = SMatrix{2,2}(Δ[1,1], Δ[2,1], Δ[1,2], Δ[2,2])
    Δs = SMatrix{2,2}(Δ[3,3], Δ[4,3], Δ[3,4], Δ[4,4])
    return (_DECOUPLED, _expm2(C, Δp, x), _expm2(C, Δs, x))
end

# p and s propagators of one interior layer, or ok = false if it couples s and p.
# Berreman's blocks for diagonal ε = (εx, εy, εz), μ = (μx, μy, μz):
#   p: [0  μy − k∥²/εz; εx  0]      s: [0  μx; εy − k∥²/μz  0]
# Rotated layers and full μ tensors go through `_tensor_blocks`.
function _interior_blocks(::Type{C}, layer, λ, k_par, k0, μ_i)::_BlockResult{C} where {C}
    ns = get_refractive_indices(layer, λ)
    _fits(C, ns...) || return (_NARROW, zero(SMatrix{2,2,C,4}), zero(SMatrix{2,2,C,4}))
    isrotated(layer) && return _tensor_blocks(C, layer, λ, k_par, k0, μ_i, ns)
    mu = _diagonal_mu(layer, λ, μ_i)
    mu === nothing && return _tensor_blocks(C, layer, λ, k_par, k0, μ_i, ns)
    k2 = k_par^2
    μx, μy, μz = mu
    nx, ny, nz = ns
    εx = dielectric_constant(nx); εy = dielectric_constant(ny); εz = dielectric_constant(nz)
    x = k0 * layer.thickness
    ap, bp = μy - k2 / εz, εx
    as, bs = μx, εy - k2 / μz
    if !isanisotropic(layer) && μx == μy == μz
        # Isotropic: both blocks share q² = με − k∥², so one cos/sin serves both.
        c, s = _cos_sinc(x, sqrt(ap * bp))
        return (_DECOUPLED, _block_from(C, ap, bp, c, s), _block_from(C, as, bs, c, s))
    end
    return (_DECOUPLED, _block_propagator(C, ap, bp, x), _block_propagator(C, as, bs, x))
end

"""
    _propagate_core_2x2([C,] λ, layers; θ=0.0, μ=1.0, sheets=nothing)

2×2 propagation core for decoupled stacks. Returns `(M_sys, S)` in the same basis
as [`_propagate_core_exp`](@ref), or `nothing` when a per-point check fails (an
interior layer whose Berreman Δ couples s and p at this λ and θ, or a sheet with
σxy ≠ 0), in which case the caller falls back to the 4×4 path. Assumes
`_decoupled_structure(layers)`. `C` is the complex element type of the
calculation (see [`_stack_ctype`](@ref)); it is computed when omitted. With a `C`
too narrow for some refractive index, returns `_WIDEN`.

Per interior layer with diagonal ε = (εx, εy, εz) and μ = (μx, μy, μz), Berreman's
blocks are p: [0  μy − k∥²/εz; εx  0] and s: [0  μx; εy − k∥²/μz  0]. A rotated
layer or full μ tensor uses the two diagonal blocks of its Δ, exponentiated in
closed form.
"""
_propagate_core_2x2(λ, layers; θ=0.0, μ=1.0, sheets=nothing) =
    _propagate_core_2x2(_stack_ctype(λ, θ, μ, layers, sheets), λ, layers; θ=θ, μ=μ, sheets=sheets)

function _propagate_core_2x2(::Type{C}, λ, layers; θ=0.0, μ=1.0, sheets=nothing) where {C}
    N = length(layers)
    ok1, ε_1 = _boundary_epsilon(C, layers[1], λ)
    okN, ε_N = _boundary_epsilon(C, layers[N], λ)
    (ok1 && okN) || return _WIDEN
    k_par = √ε_1 * sin(θ)
    res = _block_matrices_2x2(C, λ, layers, 1, N, k_par, μ, sheets, ε_1, ε_N)
    res isa Tuple || return res
    Mp, Ms, med_a, med_b = res
    z = zero(C)
    M_sys = @SMatrix [Mp[1,1] Mp[1,2] z z;
                      Mp[2,1] Mp[2,2] z z;
                      z z Ms[1,1] Ms[1,2];
                      z z Ms[2,1] Ms[2,2]]
    # calculate_tr(M_sys) with the zero blocks dropped: d = Mp₁₁·Ms₁₁
    tpp = 1 / Mp[1,1]; tss = 1 / Ms[1,1]
    r = SVector(Mp[2,1] * tpp, z, Ms[2,1] * tss, z)
    t = SVector(tpp, z, z, tss)
    return M_sys, _poynting_decoupled(med_a, med_b, t, r)
end

# Reduced in-plane wavevector k∥ = n_ambient·sin θ, conserved through the stack.


"""
    _block_matrices_2x2(C, λ, layers, a, b, k_par, μ, sheets)

p and s 2×2 transfer matrices (in the mode-amplitude basis, `[E⁺; E⁻]_a = M [E⁺; E⁻]_b`)
of the sub-stack `layers[a:b]`, whose end layers `a` and `b` are isotropic,
non-magnetic and treated as semi-infinite, together with the two boundary media
(`_IsoMedium`). `ε_a`, `ε_b` are the end layers' permittivities. Sheets are
keyed by the global interface index (sheet `i` sits between layers `i` and `i+1`).
Returns `nothing` if some layer or sheet in the range couples s and p, and `_WIDEN`
if `C` cannot hold a layer's refractive index.
"""
function _block_matrices_2x2(::Type{C}, λ, layers, a, b, k_par, μ, sheets, ε_a, ε_b) where {C}
    k0 = (2π * c_0 / λ) / c_0
    no_sheets = sheets === nothing || isempty(sheets)

    Bp = SMatrix{2,2,C}(I)
    Bs = SMatrix{2,2,C}(I)
    if !no_sheets && haskey(sheets, a)
        blocks = _sheet_blocks(C, sheets[a], λ)
        blocks === nothing && return nothing
        Bp = Bp * blocks[1]; Bs = Bs * blocks[2]
    end
    for i in a+1:b-1
        # Function barrier: a stack mixing dispersion types is a Vector of abstract
        # Layer, so dispatch once per layer rather than once per field access.
        status, P, Sb = _interior_blocks(C, layers[i], λ, k_par, k0, μ)::_BlockResult{C}
        status === _DECOUPLED || return status === _NARROW ? _WIDEN : nothing
        Bp = Bp * P
        Bs = Bs * Sb
        if !no_sheets && haskey(sheets, i)
            blocks = _sheet_blocks(C, sheets[i], λ)
            blocks === nothing && return nothing
            Bp = Bp * blocks[1]; Bs = Bs * blocks[2]
        end
    end

    med_a = _iso_medium(C, ε_a, k_par, μ)
    med_b = _iso_medium(C, ε_b, k_par, μ)

    # The isotropic dynamical matrices are block-diagonal too (p modes 1, 3 carry no
    # Eᵧ or Hₓ; s modes 2, 4 no Eₓ, E_z or Hᵧ), so D₁⁻¹·B·D_N splits into two 2×2
    # products. In the swapped basis of `_propagate_core_exp` the p block lands in
    # rows/columns 1:2 and the s block in 3:4.
    Mp = inv(med_a.Dp) * Bp * med_b.Dp
    Ms = inv(med_a.Ds) * Bs * med_b.Ds
    return Mp, Ms, med_a, med_b
end

const _BoundaryModes{C} = Tuple{SMatrix{4,4,C,16},SMatrix{4,3,C,12},SVector{4,C}}

# Permittivity of an isotropic boundary layer: the only layer-dependent (hence, for a
# Vector of abstract Layer, dynamically dispatched) step of the boundary modes, kept
# small so the dispatch boxes one number rather than the mode matrices.
# Returns (fits, ε); see `_fits`.
# (The assertion at the call site matters: `layer` is often abstractly typed, so the
# callee's return annotation alone does not reach the caller's inference.)
@inline _boundary_epsilon(::Type{C}, layer, λ) where {C} = _boundary_epsilon_impl(C, layer, λ)::Tuple{Bool,C}
function _boundary_epsilon_impl(::Type{C}, layer, λ)::Tuple{Bool,C} where {C}
    n = first(get_refractive_indices(layer, λ))
    _fits(C, n) || return (false, zero(C))
    return (true, C(dielectric_constant(n)))
end

function _boundary_modes(::Type{C}, ε, k_par, μ)::_BoundaryModes{C} where {C}
    D, E, q = _isotropic_modes(ε, k_par, μ)
    return SMatrix{4,4,C}(D), SMatrix{4,3,C}(E), SVector{4,C}(q)
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
