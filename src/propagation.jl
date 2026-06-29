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


# Lightweight path: only computes Γ and S without allocating per-layer
# Ds, Ps, γs vectors. Used by `transfer` in tight spectral loops.
function _propagate_core(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

    N = length(layers)
    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ε_0in = dielectric_constant(nx_in)
    k_par = √(ε_0in) * sin(θ)  # reduced in-plane wavevector, = n·sinθ

    Λ_1324 = @SMatrix [1 0 0 0;
                       0 0 1 0;
                       0 1 0 0;
                       0 0 0 1]

    no_sheets = sheets === nothing || isempty(sheets)

    D_prev, _, γ_first, q_first = layer_matrices(layers[1], λ, k_par, μ)
    γ_last = γ_first
    q_last = q_first

    Γ = SMatrix{4,4,ComplexF64}(I)
    for i in 2:N
        layer = layers[i]
        D_cur, P_cur, γ_cur, q_cur = layer_matrices(layer, λ, k_par, μ)
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
    S = poynting(k_par, q_first, q_last, γ_first, γ_last, t, r, μ_in_mat, μ_out_mat)

    return Γ, S
end


"""
    _propagate_core_exp(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Matrix-exponential propagation core. Interior layers propagate via
[`layer_transfer_exp`](@ref) (no eigenmode sorting); the semi-infinite ambient and
substrate keep the eigenmode treatment in [`layer_matrices`](@ref), which is needed
for the r/t coefficients and the Poynting transmittance. Conductive sheets are
injected at their respective interfaces, mirroring `_propagate_core`.
Returns `(Γ, S)` like `_propagate_core`.
"""
function _propagate_core_exp(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

    N = length(layers)
    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ε_0in = dielectric_constant(nx_in)
    k_par = √(ε_0in) * sin(θ)
    ω = 2π * c_0 / λ

    D_1, _, γ_first, q_first = layer_matrices(layers[1], λ, k_par, μ)
    D_N, _, γ_last,  q_last  = layer_matrices(layers[N], λ, k_par, μ)

    # Interior product in the dynamical-matrix basis. Tangential fields are
    # continuous across plain interfaces, so interior layers chain directly; a
    # conductive sheet at interface (i, i+1) is injected as sheet_matrix.
    no_sheets = sheets === nothing || isempty(sheets)
    core = SMatrix{4,4,ComplexF64}(I)
    if !no_sheets && haskey(sheets, 1)
        core = core * sheet_matrix(sheets[1], λ)
    end
    for i in 2:N-1
        core = core * layer_transfer_exp(layers[i], λ, k_par, ω, μ)
        if !no_sheets && haskey(sheets, i)
            core = core * sheet_matrix(sheets[i], λ)
        end
    end
    core = core * D_N

    Γ = (_Λ1324 \ (D_1 \ core)) * _Λ1324
    r, R, t, T = calculate_tr(Γ)
    μ_in_mat  = ismagnetic(layers[1])   ? get_permeability(layers[1],   λ) : SMatrix{3,3,ComplexF64}(μ*I)
    μ_out_mat = ismagnetic(layers[end]) ? get_permeability(layers[end], λ) : SMatrix{3,3,ComplexF64}(μ*I)
    S = poynting(k_par, q_first, q_last, γ_first, γ_last, t, r, μ_in_mat, μ_out_mat)

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
    k_par = √(ε_0in) * sin(θ)

    Λ_1324 = @SMatrix [1 0 0 0;
                       0 0 1 0;
                       0 1 0 0;
                       0 0 0 1]

    no_sheets = sheets === nothing || isempty(sheets)

    D_1, P_1, γ_1, q_1 = layer_matrices(layers[1], λ, k_par, μ)
    Ds = Vector{typeof(D_1)}(undef, N)
    Ps = Vector{typeof(P_1)}(undef, N)
    γs = Vector{typeof(γ_1)}(undef, N)
    qs = Vector{typeof(q_1)}(undef, N)
    Ds[1] = D_1; Ps[1] = P_1; γs[1] = γ_1; qs[1] = q_1

    Γ = SMatrix{4,4,ComplexF64}(I)
    for i in 2:N
        D_i, P_i, γ_i, q_i = layer_matrices(layers[i], λ, k_par, μ)
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
    S = poynting(k_par, q_1, qs[N], γ_1, γs[N], t, r, μ_in_mat, μ_out_mat)

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
