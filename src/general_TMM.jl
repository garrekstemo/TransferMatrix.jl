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
injected at their respective interfaces, mirroring `_propagate_core`.
Returns `(Γ, S)` like `_propagate_core`.
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

# Numerical backend

`method` selects how interior layers are propagated:

- `:exp` (default) computes each interior layer's transfer matrix as the matrix
  exponential of the Berreman Δ matrix (see [`layer_transfer_exp`](@ref)). It needs
  no eigenmode sorting and is degeneracy-immune, so it handles near-degenerate and
  mixed propagating/evanescent interior layers that `:eig` cannot.
- `:eig` is the eigenmode/dynamical-matrix path, retained as a cross-check.

Both agree to ~1e-12 on all supported cases. The semi-infinite ambient and substrate
use the eigenmode path in either mode, so the anisotropic-ambient (#71) and
anisotropic-substrate (#107) boundary limitations are unaffected by `method`.

Mackay & Lakhtakia, 2020, https://doi.org/10.1007/978-3-031-02022-3 ;
Higham, 2005, https://doi.org/10.1137/04061101X

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
function transfer(λ, layers; θ=0.0, μ=1.0, sheets=nothing, validate::Bool=false, basis::Symbol=:linear, method::Symbol=:exp)
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


function _sweep_spectra(outer_vals, inner_vals, ::Val{B}; threads::Bool=true, verbose::Bool=false, make_layers, angle_for, sheets=nothing, method::Symbol=:eig) where {B}
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
            result = transfer(inner_vals[j], layers_i; θ=θ, sheets=sheets, basis=B, method=method)
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
- `method`: propagation backend, `:exp` (default) or `:eig` (see [`transfer`](@ref))

# Units
- Wavelengths: μm (micrometers) recommended
- Angles: radians
"""
function sweep_angle(λs, θs, layers; sheets=nothing, threads::Bool=true, verbose::Bool=false, basis::Symbol=:linear, method::Symbol=:exp)
    λs = _to_wavelength_um.(λs)
    θs = _to_radians.(θs)
    sd = sheets === nothing ? nothing : _sheets_dict(sheets)
    _validate_sheet_indices(sd, length(layers))
    return _sweep_spectra(θs, λs, Val(basis); threads=threads, verbose=verbose,
        make_layers = _ -> layers,
        angle_for = i -> θs[i],
        sheets = sd, method = method)
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
- `method`: propagation backend, `:exp` (default) or `:eig` (see [`transfer`](@ref))

# Units
- Wavelengths and thicknesses: μm (micrometers) recommended
- Angle: radians
"""
function sweep_thickness(λs, ts, layers, t_index::Int; θ=0.0, sheets=nothing, threads::Bool=true, verbose::Bool=false, basis::Symbol=:linear, method::Symbol=:exp)
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
        sheets = sd, method = method)
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
