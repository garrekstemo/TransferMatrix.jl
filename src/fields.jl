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
