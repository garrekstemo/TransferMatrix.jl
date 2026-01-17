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

struct ElectricField
    z::Vector{Float64}
    p::Matrix{ComplexF64}
    s::Matrix{ComplexF64}
    boundaries::Vector{Float64}

    function ElectricField(z, p, s, boundaries)
        size(p, 2) == length(z) || throw(ArgumentError("p field columns must match z length"))
        size(s, 2) == length(z) || throw(ArgumentError("s field columns must match z length"))
        new(z, p, s, boundaries)
    end
end

struct Spectra
    Rpp::Matrix{Float64}
    Rss::Matrix{Float64}
    Tpp::Matrix{Float64}
    Tss::Matrix{Float64}
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
    poynting(ξ, q_in, q_out, γ_in, γ_out, t_coefs, r_coefs)

Calculate the Poynting vector from wavevectors ``q``,
componments of the electric field γ, and transmission
and reflection coefficients.
"""
function poynting(ξ, q_in, q_out, γ_in, γ_out, t_coefs, r_coefs)

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

    E_out_p = t_coefs[1] * γ_out[1, :] + t_coefs[2] * γ_out[2, :]
    E_out_s = t_coefs[3] * γ_out[1, :] + t_coefs[4] * γ_out[2, :]

    E_ref_p = r_coefs[1] * γ_in[3, :] + r_coefs[2] * γ_in[4, :]
    E_ref_s = r_coefs[3] * γ_in[3, :] + r_coefs[4] * γ_in[4, :]

    S_in_p = real(0.5 * E_forward_in_p × conj(k_in[1, :] × E_forward_in_p))
    S_in_s = real(0.5 * E_forward_in_s × conj(k_in[2, :] × E_forward_in_s))

    k_out = @MMatrix zeros(ComplexF64, 4, 3)
    k_out[:, 1] .= ξ

    for (i, q_i) in enumerate(q_out)
        k_out[i, 3] = q_i
    end

    k_out ./= c_0
    k_out = SMatrix(k_out)

    S_out_p  = real(0.5 * E_out_p × conj(k_out[1, :] × E_out_p))
    S_out_s  = real(0.5 * E_out_s × conj(k_out[2, :] × E_out_s))
    S_refl_p = real(0.5 * E_ref_p × conj(k_out[3, :] × E_ref_p))
    S_refl_s = real(0.5 * E_ref_s × conj(k_out[4, :] × E_ref_s))

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
function evaluate_birefringence(Ψ, S, t_modes, r_modes)

    C_q1 = abs_ratio(S[1, t_modes[1]], S[2, t_modes[1]])
    C_q2 = abs_ratio(S[1, t_modes[2]], S[2, t_modes[2]])
    
    # Note: isapprox(NaN, NaN) = false,
    # which is important in the case that both Sx and Sy are zero.
    if isapprox(C_q1, C_q2)

        if C_q2 > C_q1
            reverse!(t_modes)
        end

        C_q3 = abs_ratio(S[1, r_modes[1]], S[2, r_modes[1]])
        C_q4 = abs_ratio(S[1, r_modes[2]], S[2, r_modes[2]])

        if C_q4 > C_q3
            reverse!(r_modes)
        end
    else
        C_q1 = abs_ratio(Ψ[1, t_modes[1]], Ψ[3, t_modes[1]])
        C_q2 = abs_ratio(Ψ[1, t_modes[2]], Ψ[3, t_modes[2]])

        if C_q2 > C_q1
            reverse!(t_modes)
        end

        C_q3 = abs_ratio(Ψ[1, r_modes[1]], Ψ[3, r_modes[1]])
        C_q4 = abs_ratio(Ψ[1, r_modes[2]], Ψ[3, r_modes[2]])

        if C_q4 > C_q3
            reverse!(r_modes)
        end
    end

    return t_modes, r_modes
end


"""
Ratio of the absolution squares of two components
used to evaluate if a material is birefringent.
"""
abs_ratio(a, b) = abs2(a) / (abs2(a) + abs2(b))


"""
    propagate(λ, layers; θ=0.0, μ=1.0)

Calculate the transfer matrix for the entire structure,
as well as the Poynting vector for the structure.

# Arguments
- `λ`: Wavelength
- `layers`: Vector of `Layer` objects representing the stack
- `θ`: Angle of incidence in radians (default: 0.0, normal incidence)
- `μ`: Relative magnetic permeability (default: 1.0, non-magnetic)
"""
function propagate(λ, layers; θ=0.0, μ=1.0)

    first_layer = layers[1]
    last_layer = layers[end]

    n_in = first_layer.dispersion(λ)
    ε_0in = dielectric_constant(n_in)
    ξ = √(ε_0in) * sin(θ)

    Λ_1324 = @SMatrix [1 0 0 0;
                       0 0 1 0;
                       0 1 0 0;
                       0 0 0 1]
   
    D_0, P_0, γ_0, q_0 = layer_matrices(first_layer, λ, ξ, μ)
    D_f, P_f, γ_f, q_f = layer_matrices(last_layer, λ, ξ, μ)

    # Preallocate arrays with known size
    n_layers = length(layers)
    Ds = Vector{typeof(D_0)}(undef, n_layers)
    Ps = Vector{typeof(P_0)}(undef, n_layers)
    γs = Vector{typeof(γ_0)}(undef, n_layers)

    Ds[1] = D_0
    Ps[1] = P_0
    γs[1] = γ_0

    Γ = SMatrix{4,4,ComplexF64}(I)
    for (i, layer) in enumerate(layers[2:end - 1])
        D_i, P_i, γ_i, q_i = layer_matrices(layer, λ, ξ, μ)
        # Prefer solves over inv() for stability and fewer allocations.
        T_i = D_i * (P_i(layer.thickness) / D_i)
        Γ *= T_i
        Ds[i + 1] = D_i
        Ps[i + 1] = P_i
        γs[i + 1] = γ_i
    end

    Ds[n_layers] = D_f
    Ps[n_layers] = P_f
    γs[n_layers] = γ_f

    Γ = (Λ_1324 \ (D_0 \ (Γ * D_f))) * Λ_1324
    r, R, t, T = calculate_tr(Γ)
    S = poynting(ξ, q_0, q_f, γ_0, γ_f, t, r)

    return Γ, S, Ds, Ps, γs
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
    calculate_tr(λ, layers; θ=0.0, μ=1.0, validate=false)

Calculate the transmittance and reflectance of a layered structure.

Returns `(Tpp, Tss, Rpp, Rss)` where:
- `Tpp`, `Tss`: p- and s-polarized transmittance
- `Rpp`, `Rss`: p- and s-polarized reflectance

Transmittance is calculated via Poynting vectors for accurate energy flow.
Reflectance is calculated directly from transfer matrix elements.

# Arguments
- `λ`: Wavelength in μm (must match units used for layer thicknesses)
- `layers`: Vector of `Layer` objects representing the stack
- `θ`: Angle of incidence in radians (default: 0.0, normal incidence)
- `μ`: Relative magnetic permeability (default: 1.0, non-magnetic)
- `validate`: Check energy conservation R + T ≈ 1 for non-absorbing media (default: false)

# Wave Propagation Convention
- Light propagates in the **+z direction** (from first layer toward last layer)
- The first and last layers are treated as semi-infinite media
- θ is measured from the surface normal (z-axis)

# Units
- Wavelength and thicknesses: μm (micrometers) recommended
- Angle: radians
- Transmittance/Reflectance: dimensionless (0 to 1)

# Energy Conservation Validation
When `validate=true`, the function checks if all layers are non-absorbing (imaginary part
of refractive index < 1e-10) and verifies that R + T ≈ 1 (within tolerance of 1e-6).
A warning is issued if energy conservation is violated.
"""
function calculate_tr(λ, layers; θ=0.0, μ=1.0, validate::Bool=false)

    Γ, S, Ds, Ps, γs = propagate(λ, layers; θ=θ, μ=μ)
    r, R, t, T = calculate_tr(Γ)
    Tpp, Tss, Rpp_, Rss_ = calculate_tr(S)
    Rpp = R[1]
    Rss = R[2]

    if validate
        _validate_energy_conservation(λ, layers, Tpp, Tss, Rpp, Rss)
    end

    return Tpp, Tss, Rpp, Rss
end


"""
    _validate_energy_conservation(λ, layers, Tpp, Tss, Rpp, Rss; atol=1e-6, k_threshold=1e-10)

Check energy conservation (R + T ≈ 1) for non-absorbing media.
Issues a warning if energy is not conserved within tolerance.

Internal function called by `calculate_tr` when `validate=true`.
"""
function _validate_energy_conservation(λ, layers, Tpp, Tss, Rpp, Rss; atol=1e-6, k_threshold=1e-10)
    # Check if all layers are non-absorbing
    is_lossless = all(layers) do layer
        n = layer.dispersion(λ)
        abs(imag(n)) < k_threshold
    end

    if is_lossless
        sum_p = Tpp + Rpp
        sum_s = Tss + Rss

        if !isapprox(sum_p, 1.0; atol=atol)
            @warn "Energy conservation violated for p-polarization" Tpp Rpp sum=sum_p expected=1.0 deviation=abs(sum_p - 1.0)
        end

        if !isapprox(sum_s, 1.0; atol=atol)
            @warn "Energy conservation violated for s-polarization" Tss Rss sum=sum_s expected=1.0 deviation=abs(sum_s - 1.0)
        end
    end

    return nothing
end


"""
    sweep_angle(λs, θs, layers; threads=true, verbose=false)

Calculate transmittance/reflectance spectra over wavelength and angle of incidence.

Returns a `Spectra` struct with fields `Rpp`, `Rss`, `Tpp`, `Tss`, each a matrix
of size `(length(θs), length(λs))`.

# Arguments
- `λs`: Vector of wavelengths in μm
- `θs`: Vector of angles of incidence in radians
- `layers`: `AbstractVector{<:Layer}` representing the stack
- `threads`: Enable multithreading (default: true)
- `verbose`: Print thread count info (default: false)

# Units
- Wavelengths: μm (micrometers) recommended
- Angles: radians
"""
function _sweep_spectra(outer_vals, inner_vals; threads::Bool=true, verbose::Bool=false, make_layers, angle_for)
    Tpp = Array{Float64}(undef, length(outer_vals), length(inner_vals))
    Tss = Array{Float64}(undef, length(outer_vals), length(inner_vals))
    Rpp = Array{Float64}(undef, length(outer_vals), length(inner_vals))
    Rss = Array{Float64}(undef, length(outer_vals), length(inner_vals))

    if verbose
        println("Threads: ", Threads.nthreads())
    end

    function compute_row(i)
        layers_i = make_layers(i)
        θ = angle_for(i)
        for j in eachindex(inner_vals)
            Tpp_, Tss_, Rpp_, Rss_ = calculate_tr(inner_vals[j], layers_i; θ=θ)
            Tpp[i, j] = Tpp_
            Tss[i, j] = Tss_
            Rpp[i, j] = Rpp_
            Rss[i, j] = Rss_
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

    return Spectra(Rpp, Rss, Tpp, Tss)
end

function sweep_angle(λs, θs, layers; threads::Bool=true, verbose::Bool=false)
    return _sweep_spectra(θs, λs; threads=threads, verbose=verbose,
        make_layers = _ -> layers,
        angle_for = i -> θs[i])
end


"""
    sweep_thickness(λs, ts, layers, t_index; θ=0.0, threads=true, verbose=false)

Sweep the thickness of a specific layer and calculate transmittance/reflectance spectra.

Returns a `Spectra` struct with fields `Rpp`, `Rss`, `Tpp`, `Tss`, each a matrix
of size `(length(ts), length(λs))`.

# Arguments
- `λs`: Vector of wavelengths in μm
- `ts`: Vector of thicknesses in μm to sweep
- `layers`: `AbstractVector{<:Layer}` representing the stack
- `t_index`: Index of the layer (1-based) whose thickness to vary
- `θ`: Angle of incidence in radians (default: 0.0, normal incidence)
- `threads`: Enable multithreading (default: true)
- `verbose`: Print thread count info (default: false)

# Units
- Wavelengths and thicknesses: μm (micrometers) recommended
- Angle: radians
"""
function sweep_thickness(λs, ts, layers, t_index::Int; θ=0.0, threads::Bool=true, verbose::Bool=false)
    # Pre-allocate mutable layer vector to avoid repeated slicing/concatenation
    layers_mut = collect(layers)
    dispersion_func = layers[t_index].dispersion

    return _sweep_spectra(ts, λs; threads=threads, verbose=verbose,
        make_layers = i -> begin
            layers_mut[t_index] = Layer(dispersion_func, ts[i])
            layers_mut
        end,
        angle_for = _ -> θ)
end

@deprecate angle_resolved(λs, θs, layers; kwargs...) sweep_angle(λs, θs, layers; kwargs...)
@deprecate tune_thickness(λs, ts, layers, t_index::Int, θ=0.0; kwargs...) sweep_thickness(λs, ts, layers, t_index; θ=θ, kwargs...)


"""
    electric_field(λ, layers; θ=0.0, μ=1.0, dz=0.001)

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

# Wave Propagation Convention
- Light propagates in the **+z direction** (from first layer toward last layer)
- z = 0 is at the first interface (between layer 1 and layer 2)
- Negative z values are inside the first (incident) layer
- θ is measured from the surface normal (z-axis)

# Units
- All lengths (λ, thickness, dz, z): μm (micrometers) recommended
- Angle: radians
- Electric field: arbitrary units (normalized to incident field)
"""
function electric_field(λ, layers; θ=0.0, μ=1.0, dz=0.001)

    Γ, S, Ds, Ps, γs = propagate(λ, layers; θ=θ, μ=μ)
    r, R, t, T = calculate_tr(Γ)
    first_layer = layers[1]
    last_layer = layers[end]

    Eplus_p = zeros(ComplexF64, length(layers), 4)
    Eminus_p = zeros(ComplexF64, length(layers), 4)

    Eplus_s = zeros(ComplexF64, length(layers), 4)
    Eminus_s = zeros(ComplexF64, length(layers), 4)

    Eplus_p[end, :] = [t[1], t[2], 0, 0]
    Eplus_s[end, :] = [t[3], t[4], 0, 0]
    
    P_f = Ps[end]
    # Avoid inv() here too; P_f is diagonal so the solve is cheap and stable.
    Eminus_p[end, :] = P_f(last_layer.thickness) \ Eplus_p[end, :]
    Eminus_s[end, :] = P_f(last_layer.thickness) \ Eplus_s[end, :]

    D_i = Ds[end]

    for l in reverse(eachindex(layers))
        if l >= 2
            layer = layers[l - 1]
            D_prev = Ds[l - 1]
            P_prev = Ps[l - 1]
            L_i = Ds[l - 1] \ D_i

            Eminus_p[l - 1, :] = L_i * Eplus_p[l, :]
            Eminus_s[l - 1, :] = L_i * Eplus_s[l, :]
            Eplus_p[l - 1, :] = P_prev(layer.thickness) * Eminus_p[l - 1, :]
            Eplus_s[l - 1, :] = P_prev(layer.thickness) * Eminus_s[l - 1, :]
            
            D_i = D_prev
        end
    end

    interface_positions, total_thickness = find_bounds(layers)
    interface_positions .-= first_layer.thickness

    zs = range(-first_layer.thickness, interface_positions[end], step=dz)

    field = zeros(ComplexF64, 6, length(zs))

    i = 1
    for (j, z) in enumerate(zs)
        if i < length(layers) && z > interface_positions[i]
            i += 1
        end

        P_i = Ps[i]
        field_p = P_i(-(z - interface_positions[i])) * Eminus_p[i, :]
        field_s = P_i(-(z - interface_positions[i])) * Eminus_s[i, :]

        @views field[1:3, j] = field_p[1] * γs[i][1, :] + field_p[2] * γs[i][2, :] + field_p[3] * γs[i][3, :] + field_p[4] * γs[i][4, :]
        @views field[4:6, j] = field_s[1] * γs[i][1, :] + field_s[2] * γs[i][2, :] + field_s[3] * γs[i][3, :] + field_s[4] * γs[i][4, :]
    end

    return ElectricField(zs, field[1:3, :], field[4:6, :], interface_positions[1:end - 1])
end
