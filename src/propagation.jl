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


# Lightweight path: only computes M_sys and S without allocating per-layer
# Ds, Ps, E_modes_per_layer vectors. Used by `transfer` in tight spectral loops.
function _propagate_core(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

    N = length(layers)
    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ε_0in = dielectric_constant(nx_in)
    k_par = √(ε_0in) * sin(θ)  # reduced in-plane wavevector, = n·sinθ

    swap23 = @SMatrix [1 0 0 0;
                       0 0 1 0;
                       0 1 0 0;
                       0 0 0 1]

    no_sheets = sheets === nothing || isempty(sheets)

    D_prev, _, E_modes_first, q_first = layer_matrices(layers[1], λ, k_par, μ)
    E_modes_last = E_modes_first
    q_last = q_first

    M_sys = SMatrix{4,4,ComplexF64}(I)
    for i in 2:N
        layer = layers[i]
        D_cur, P_cur, E_modes_cur, q_cur = layer_matrices(layer, λ, k_par, μ)
        if no_sheets || !haskey(sheets, i - 1)
            L = D_prev \ D_cur                                  # interface (i-1, i)
        else
            L = D_prev \ (sheet_matrix(sheets[i - 1], λ) * D_cur)
        end
        M_sys *= L                                              # first ⇒ D₀⁻¹D₂ ; last ⇒ D_{N-1}⁻¹D_f
        if i < N
            M_sys *= P_cur(layer.thickness)                     # propagate interior layer i
        end
        D_prev = D_cur
        if i == N
            E_modes_last = E_modes_cur
            q_last = q_cur
        end
    end

    M_sys = (swap23 \ M_sys) * swap23
    r, R, t, T = calculate_tr(M_sys)
    μ_in_mat  = ismagnetic(layers[1])   ? get_permeability(layers[1],   λ) : SMatrix{3,3,ComplexF64}(μ*I)
    μ_out_mat = ismagnetic(layers[end]) ? get_permeability(layers[end], λ) : SMatrix{3,3,ComplexF64}(μ*I)
    S = poynting(k_par, q_first, q_last, E_modes_first, E_modes_last, t, r, μ_in_mat, μ_out_mat)

    return M_sys, S
end


"""
    _propagate_core_exp(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Matrix-exponential propagation core. Interior layers propagate via
[`layer_transfer_exp`](@ref) (no eigenmode sorting); the semi-infinite ambient and
substrate keep the eigenmode treatment in [`layer_matrices`](@ref), which is needed
for the r/t coefficients and the Poynting transmittance. Conductive sheets are
injected at their respective interfaces, mirroring `_propagate_core`.
Returns `(M_sys, S)` like `_propagate_core`.
"""
_propagate_core_exp(λ, layers; θ=0.0, μ=1.0, sheets=nothing) =
    _propagate_core_exp(_stack_ctype(λ, θ, μ, layers, sheets), λ, layers; θ=θ, μ=μ, sheets=sheets)

function _propagate_core_exp(::Type{C}, λ, layers; θ=0.0, μ=1.0, sheets=nothing) where {C}

    N = length(layers)
    ok, ε_1 = _boundary_epsilon(C, layers[1], λ)
    ok || return _WIDEN
    k_par = √ε_1 * sin(θ)
    ω = 2π * c_0 / λ

    ok1, D_1, E_modes_first, q_first = _semi_infinite_modes(C, layers[1], λ, k_par, μ)::_SemiInfinite{C}
    okN, D_N, E_modes_last,  q_last  = _semi_infinite_modes(C, layers[N], λ, k_par, μ)::_SemiInfinite{C}
    (ok1 && okN) || return _WIDEN

    # Interior product in the dynamical-matrix basis. Tangential fields are
    # continuous across plain interfaces, so interior layers chain directly; a
    # conductive sheet at interface (i, i+1) is injected as sheet_matrix.
    no_sheets = sheets === nothing || isempty(sheets)
    core = SMatrix{4,4,C}(I)
    if !no_sheets && haskey(sheets, 1)
        core = core * sheet_matrix(sheets[1], λ)
    end
    for i in 2:N-1
        ok, T_i = _layer_transfer_exp(C, layers[i], λ, k_par, ω, μ)::Tuple{Bool,SMatrix{4,4,C,16}}
        ok || return _WIDEN
        core = core * T_i
        if !no_sheets && haskey(sheets, i)
            core = core * sheet_matrix(sheets[i], λ)
        end
    end
    core = core * D_N

    M_sys = (_swap23 \ (D_1 \ core)) * _swap23
    r, R, t, T = calculate_tr(M_sys)
    μ_in_mat  = _ambient_mu(C, layers[1], λ, μ)::SMatrix{3,3,C,9}
    μ_out_mat = _ambient_mu(C, layers[N], λ, μ)::SMatrix{3,3,C,9}
    S = poynting(k_par, q_first, q_last, E_modes_first, E_modes_last, t, r, μ_in_mat, μ_out_mat)

    return M_sys, S
end

# Function barriers with concrete return types for the 4×4 path, so a Vector of
# abstract Layer costs one dispatch per layer and the product stays type-stable.
# Each returns `fits` first (see `_fits`).
function _layer_transfer_exp(::Type{C}, layer, λ, k_par, ω, μ)::Tuple{Bool,SMatrix{4,4,C,16}} where {C}
    ns = get_refractive_indices(layer, λ)
    _fits(C, ns...) || return (false, zero(SMatrix{4,4,C,16}))
    return (true, SMatrix{4,4,C}(layer_transfer_exp(layer, λ, k_par, ω, μ, ns)))
end

_ambient_mu(::Type{C}, layer, λ, μ) where {C} =
    (ismagnetic(layer) ? SMatrix{3,3,C}(get_permeability(layer, λ)) : SMatrix{3,3,C}(μ * I))::SMatrix{3,3,C,9}

# Dynamical matrix, E modes and wavevectors of a semi-infinite ambient or substrate.
# Isotropic, non-magnetic media use the closed-form modes of the 2×2 path (the same
# basis as the eigen-solve, without it, so this also works with dual numbers);
# anything else is diagonalized by `layer_matrices`.
const _SemiInfinite{C} = Tuple{Bool,SMatrix{4,4,C,16},SMatrix{4,3,C,12},SVector{4,C}}

function _semi_infinite_modes(::Type{C}, layer, λ, k_par, μ)::_SemiInfinite{C} where {C}
    if _isotropic_boundary(layer)
        ok, ε = _boundary_epsilon(C, layer, λ)
        return (ok, _boundary_modes(C, ε, k_par, μ)...)
    end
    _fits(C, get_refractive_indices(layer, λ)...) ||
        return (false, zero(SMatrix{4,4,C,16}), zero(SMatrix{4,3,C,12}), zero(SVector{4,C}))
    D, _, E, q = layer_matrices(layer, λ, k_par, μ)
    return (true, SMatrix{4,4,C}(D), SMatrix{4,3,C}(E), SVector{4,C}(q))
end

# Dispatch between the cores: eigenmode (:eig), matrix exponential (:exp), the 2×2
# decoupled path (:decoupled, see decoupled.jl), and :auto, which takes the 2×2
# path whenever s and p decouple at this point and the 4×4 :exp path otherwise.
_propagate(::Val{:exp}, λ, layers; kwargs...) = _propagate_core_exp(λ, layers; kwargs...)
_propagate(::Val{:eig}, λ, layers; kwargs...) = _propagate_core(λ, layers; kwargs...)
_propagate(m::Val, λ, layers; θ=0.0, μ=1.0, sheets=nothing) =
    _propagate(m, _stack_ctype(λ, θ, μ, layers, sheets), λ, layers; θ=θ, μ=μ, sheets=sheets)

# Typed methods: `C` is the complex element type of the calculation (`_stack_ctype`).
# The eigenmode path diagonalizes every layer with `eigen`, which has no dual-number
# method, so it ignores `C` and runs in ComplexF64.
_propagate(::Val{:exp}, ::Type{C}, λ, layers; kwargs...) where {C} = _propagate_core_exp(C, λ, layers; kwargs...)
_propagate(::Val{:eig}, ::Type{C}, λ, layers; kwargs...) where {C} = _propagate_core(λ, layers; kwargs...)
function _propagate(::Val{:auto}, ::Type{C}, λ, layers; kwargs...) where {C}
    if _decoupled_structure(layers)
        res = _propagate_core_2x2(C, λ, layers; kwargs...)
        res === nothing || return res        # (M_sys, S) or _WIDEN
    end
    return _propagate_core_exp(C, λ, layers; kwargs...)
end
function _propagate(::Val{:decoupled}, ::Type{C}, λ, layers; kwargs...) where {C}
    res = _decoupled_structure(layers) ? _propagate_core_2x2(C, λ, layers; kwargs...) : nothing
    res === nothing && throw(ArgumentError(
        "method=:decoupled needs a stack in which s and p decouple (isotropic, " *
        "non-magnetic ambient and substrate; no interior layer coupling them at this λ " *
        "and θ, i.e. y a principal axis of every ε and μ; " *
        "sheets with σxy = σyx = 0). Use method=:auto to fall back to the 4×4 path."))
    return res
end
_propagate(::Val{M}, ::Type, λ, layers; kwargs...) where {M} =
    throw(ArgumentError("method must be :auto, :exp, :eig or :decoupled, got :$(M)"))


"""
    _propagate_full(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Internal full transfer-matrix pass. Returns `(M_sys, S, Ds, Ps, E_modes_per_layer, qs)` — like
[`propagate`](@ref) but also returns the per-layer eigenvalue vectors `qs`,
needed for magnetic-field reconstruction. Supports conductive sheets (Task: sheets).
"""
function _propagate_full(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

    N = length(layers)
    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ε_0in = dielectric_constant(nx_in)
    k_par = √(ε_0in) * sin(θ)

    swap23 = @SMatrix [1 0 0 0;
                       0 0 1 0;
                       0 1 0 0;
                       0 0 0 1]

    no_sheets = sheets === nothing || isempty(sheets)

    D_1, P_1, E_modes_1, q_1 = layer_matrices(layers[1], λ, k_par, μ)
    Ds = Vector{typeof(D_1)}(undef, N)
    Ps = Vector{typeof(P_1)}(undef, N)
    E_modes_per_layer = Vector{typeof(E_modes_1)}(undef, N)
    qs = Vector{typeof(q_1)}(undef, N)
    Ds[1] = D_1; Ps[1] = P_1; E_modes_per_layer[1] = E_modes_1; qs[1] = q_1

    M_sys = SMatrix{4,4,ComplexF64}(I)
    for i in 2:N
        D_i, P_i, E_modes_i, q_i = layer_matrices(layers[i], λ, k_par, μ)
        Ds[i] = D_i; Ps[i] = P_i; E_modes_per_layer[i] = E_modes_i; qs[i] = q_i
        if no_sheets || !haskey(sheets, i - 1)
            L = Ds[i - 1] \ D_i
        else
            L = Ds[i - 1] \ (sheet_matrix(sheets[i - 1], λ) * D_i)
        end
        M_sys *= L
        if i < N
            M_sys *= P_i(layers[i].thickness)
        end
    end

    M_sys = (swap23 \ M_sys) * swap23
    r, R, t, T = calculate_tr(M_sys)
    μ_in_mat  = ismagnetic(layers[1])   ? get_permeability(layers[1],   λ) : SMatrix{3,3,ComplexF64}(μ*I)
    μ_out_mat = ismagnetic(layers[end]) ? get_permeability(layers[end], λ) : SMatrix{3,3,ComplexF64}(μ*I)
    S = poynting(k_par, q_1, qs[N], E_modes_1, E_modes_per_layer[N], t, r, μ_in_mat, μ_out_mat)

    return M_sys, S, Ds, Ps, E_modes_per_layer, qs
end

"""
    propagate(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Calculate the transfer matrix and Poynting vector for the structure, plus the
per-layer `D`, `P`, and `E_modes` matrices used for field reconstruction. Returns the
5-tuple `(M_sys, S, Ds, Ps, E_modes_per_layer)`. See [`transfer`](@ref) for the public R/T API.
"""
propagate(λ, layers; θ=0.0, μ=1.0, sheets=nothing) =
    _propagate_full(λ, layers; θ=θ, μ=μ, sheets=sheets)[1:5]
