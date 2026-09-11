"""
    ModeAmplitudes

Per-layer Berreman mode amplitudes of a stack for unit-amplitude p- and
s-polarized incidence, as returned by [`mode_amplitudes`](@ref).

# Fields
- `p::Vector{SVector{4,ComplexF64}}`: amplitude 4-vector of each layer for p incidence
- `s::Vector{SVector{4,ComplexF64}}`: amplitude 4-vector of each layer for s incidence
- `qs::Vector{SVector{4,ComplexF64}}`: reduced normal wavevector `q = k_z/k₀` of the
  four modes of each layer (`k₀ = 2π/λ`)
- `E_modes_per_layer::Vector{SMatrix{4,3,ComplexF64}}`: unit-norm electric-field
  polarization vector of each mode (rows) in each layer
- `boundaries::Vector{Float64}`: z-positions of the `N-1` interfaces (`boundaries[1] = 0`)
- `k_par::ComplexF64`: reduced in-plane wavevector `n₁ sin θ`

# Slot convention
Slots 1 and 2 are the **forward** (+z) modes — p-like then s-like — and slots
3 and 4 the **backward** (−z) modes, again p-like then s-like. For an isotropic
layer at normal incidence only slots 1 and 3 are populated for p incidence
(`E_modes[1,:] = (1,0,0)`, `E_modes[3,:] = (−1,0,0)`) and only slots 2 and 4 for
s incidence (`E_modes[2,:] = E_modes[4,:] = (0,1,0)`); the backward p basis
vector points along −x, so a p-incidence amplitude `a₃` contributes `−a₃` to `Ex`.

# Reference plane and propagation sign
The field inside layer `i` is

```math
\\mathbf{E}_i(z) = \\sum_{m=1}^{4} a_m \\, \\mathbf{E}_m \\, e^{\\,i k_0 q_m (z - z_i)},
```

with the `exp(−iωt)` time convention (forward modes have `Re q > 0`, or
`Im q > 0` when evanescent). The reference plane `z_i` is the layer's **entry
(left) face**, `z_i = boundaries[i-1]`, for every layer `i ≥ 2` — including the
semi-infinite exit medium, whose amplitudes are `(tpp, tps, 0, 0)` for p and
`(tsp, tss, 0, 0)` for s incidence. The semi-infinite incident medium (`i = 1`)
has no entry face and is referenced to its exit face `z₁ = boundaries[1] = 0`,
where its amplitudes are `(1, 0, rpp, rps)` for p and `(0, 1, rsp, rss)` for s
incidence (see [`amplitude_coefficients`](@ref)).
"""
struct ModeAmplitudes
    p::Vector{SVector{4,ComplexF64}}
    s::Vector{SVector{4,ComplexF64}}
    qs::Vector{SVector{4,ComplexF64}}
    E_modes_per_layer::Vector{SMatrix{4,3,ComplexF64,12}}
    boundaries::Vector{Float64}
    k_par::ComplexF64
end

"""
    mode_amplitudes(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Per-layer Berreman mode amplitudes for unit-amplitude p and s incidence,
returned as a [`ModeAmplitudes`](@ref). These are the coefficients that
[`efield`](@ref) and [`hfield`](@ref) expand into field profiles; the struct
docstring states the slot ordering (forward p, forward s, backward p, backward s),
the reference plane of each layer (entry face, except the incident medium which
is referenced to `z = 0`) and the propagation sign convention
(`exp(+i k₀ q (z − z_ref))`, `exp(−iωt)` time dependence).

Arguments and units match [`efield`](@ref): `λ` in μm (or a Unitful quantity),
`θ` in radians, `sheets` as in [`transfer`](@ref). The eigenmode backend is
used, as for `efield`.

# Example
```julia
A = mode_amplitudes(1.0, layers)
A.p[1]           # (1, 0, rpp, rps) — incident medium at z = 0
A.p[end]         # (tpp, tps, 0, 0) — exit medium at the last interface
a = A.s[2]       # s incidence, layer 2, referenced to its entry face
# E_y(z) in layer 2:  a[2] exp(i k₀ q₂ (z - z₁)) + a[4] exp(i k₀ q₄ (z - z₁))
```
"""
function mode_amplitudes(λ, layers; θ=0.0, μ=1.0, sheets=nothing)
    λ = _to_wavelength_um(λ)
    θ = _to_radians(θ)
    sd = sheets === nothing ? nothing : _sheets_dict(sheets)
    _validate_sheet_indices(sd, length(layers))

    A = _mode_amplitudes(λ, layers; θ=θ, μ=μ, sd=sd)
    nlay = length(layers)
    # Layer 1 has no entry face: reference it to z = 0 (its exit face) so the
    # amplitudes read (1, 0, r, r'). Every other layer is entry-face referenced.
    p = [SVector{4,ComplexF64}(i == 1 ? A.Eminus_p[1, :] : A.Eplus_p[i, :]) for i in 1:nlay]
    s = [SVector{4,ComplexF64}(i == 1 ? A.Eminus_s[1, :] : A.Eplus_s[i, :]) for i in 1:nlay]
    qs = [SVector{4,ComplexF64}(q) for q in A.qs]
    E_modes = [SMatrix{4,3,ComplexF64,12}(E) for E in A.E_modes_per_layer]
    return ModeAmplitudes(p, s, qs, E_modes, A.interface_positions[1:end-1], ComplexF64(A.k_par))
end


"""
    LayerEnergies

Closed-form electric-energy integrals of the finite layers of a stack, as
returned by [`layer_energies`](@ref). Element `j` refers to `layers[j+1]`
(the two semi-infinite media are excluded).

# Fields
- `p::Vector{Float64}`: `∫ E† ε_H E dz` over each finite layer for p incidence
- `s::Vector{Float64}`: the same for s incidence
- `boundaries::Vector{Float64}`: z-positions of the interfaces (as in [`efield`](@ref))
"""
struct LayerEnergies
    p::Vector{Float64}
    s::Vector{Float64}
    boundaries::Vector{Float64}
end

"""
    layer_energies(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Closed-form integral of the electric energy density over every **finite**
layer of the stack, for unit-amplitude p and s incidence, returned as a
[`LayerEnergies`](@ref) (vectors of length `length(layers) - 2`; element `j`
is `layers[j+1]`).

For each finite layer with entry-face amplitudes `a` from
[`mode_amplitudes`](@ref) the field is `E(z) = Σₘ aₘ Eₘ exp(i k₀ qₘ z)` on
`0 ≤ z ≤ d`, and

```math
U = \\int_0^d \\mathbf{E}^\\dagger\\, ε_H\\, \\mathbf{E}\\, dz
  = \\sum_{m,m'} \\bar a_m a_{m'} \\left(\\mathbf{E}_m^\\dagger ε_H \\mathbf{E}_{m'}\\right)
    \\int_0^d e^{\\,i k_0 (q_{m'} - \\bar q_m) z}\\, dz ,
```

where the z-integral is evaluated analytically (it reduces to `d` when
`q_{m'} = conj(q_m)`, e.g. the self-terms of a lossless mode or the cross terms of
an evanescent pair) and `ε_H = (ε + ε†)/2` is the Hermitian part of the layer's
dielectric tensor in the lab frame — the same tensor the transfer matrix uses
(`dielectric_constant(n)` on each principal axis, rotated by the layer's Euler
angles). For a symmetric ε this is `Re ε`, so for an isotropic lossless layer
`U = n² ∫|E|² dz` and for an absorbing one `U = Re(n²) ∫|E|² dz`.

# Units
`E` is normalized to unit incident amplitude (the incident mode vector has unit
norm) and `z` is in μm, so `U` is in μm times the (dimensionless) relative
permittivity: a vacuum layer returns `∫|E|² dz` in μm. Ratios such as the
electric-energy fraction of a cavity mode, `U_gap / Σ U`, are dimensionless.

# Example
```julia
U = layer_energies(λ, layers)
f_E = U.p[gap - 1] / sum(U.p)     # gap = index of the cavity layer in `layers`
```
"""
function layer_energies(λ, layers; θ=0.0, μ=1.0, sheets=nothing)
    λ = _to_wavelength_um(λ)
    A = mode_amplitudes(λ, layers; θ=θ, μ=μ, sheets=sheets)
    k0 = 2π / λ
    N = length(layers)
    N ≥ 3 || return LayerEnergies(Float64[], Float64[], A.boundaries)
    Up = Vector{Float64}(undef, N - 2)
    Us = Vector{Float64}(undef, N - 2)
    for i in 2:N-1
        ε = SMatrix{3,3,ComplexF64}(_layer_epsilon(layers[i], λ))
        εH = (ε + ε') / 2
        d = Float64(layers[i].thickness)
        Up[i-1] = _mode_energy(A.p[i], A.E_modes_per_layer[i], A.qs[i], εH, k0, d)
        Us[i-1] = _mode_energy(A.s[i], A.E_modes_per_layer[i], A.qs[i], εH, k0, d)
    end
    return LayerEnergies(Up, Us, A.boundaries)
end

# ∫₀ᵈ exp(i k0 Δ z) dz, with the Δ → 0 limit (constant integrand) handled by a
# short series so the self-terms of lossless modes and the cross terms of an
# evanescent ± pair do not divide by zero.
function _phase_integral(Δ, k0, d)
    x = k0 * Δ * d
    if abs(x) < 1e-8
        return d * (1 + im * x / 2)
    else
        return (exp(im * x) - 1) / (im * k0 * Δ)
    end
end

# Σ_{m,m'} conj(a_m) a_m' (E_m† εH E_m') ∫ exp(i k0 (q_m' - conj(q_m)) z) dz.
# Real by construction for Hermitian εH (the double sum is a Hermitian form).
function _mode_energy(a, E_modes, q, εH, k0, d)
    U = zero(ComplexF64)
    for m in 1:4
        a[m] == 0 && continue
        Em = SVector{3,ComplexF64}(E_modes[m, 1], E_modes[m, 2], E_modes[m, 3])
        for n in 1:4
            a[n] == 0 && continue
            En = SVector{3,ComplexF64}(E_modes[n, 1], E_modes[n, 2], E_modes[n, 3])
            overlap = dot(Em, εH * En)          # E_m† εH E_n
            U += conj(a[m]) * a[n] * overlap * _phase_integral(q[n] - conj(q[m]), k0, d)
        end
    end
    return real(U)
end
