"""
    transfer(λ, layers; θ=0.0, μ=1.0, validate=false, basis=:linear, method=:exp)

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
    M_sys, S = _propagate(Val(method), λ, layers; θ=θ, μ=μ, sheets=sd)
    return _assemble(Val(basis), M_sys, S, λ, layers, sd, validate)
end

function _assemble(::Val{:linear}, M_sys, S, λ, layers, sd, validate)
    r, R, t, T = calculate_tr(M_sys)
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

function _assemble(::Val{:circular}, M_sys, S, λ, layers, sd, validate)
    r, _, t, _ = calculate_tr(M_sys)
    Tpp, Tss, _, _ = calculate_tr(S)
    return _circular_result(r, t, Tpp, Tss)
end

_assemble(::Val{B}, M_sys, S, λ, layers, sd, validate) where {B} =
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
