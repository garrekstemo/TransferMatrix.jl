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
    z::Array{Float64}
    p::Matrix{ComplexF64}
    s::Matrix{ComplexF64}
    boundaries::Array{Float64}
end

struct AngleResolvedResult
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
DOI: 10.1364/JOSA.62.000502
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

Li et al., 1988,
DOI: 10.1364/AO.27.001334

and the use of the Poynting vector is from

Passler et al., 2017,
DOI: 10.1364/JOSAB.34.002128

and

Passler et al., 2019,
DOI: 10.1364/JOSAB.36.003246
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
    propagate(λ, layers, θ, μ)

Calculate the transfer matrix for the entire structure,
as well as the Poynting vector for the structure.
"""
function propagate(λ, layers, θ, μ)

    first_layer = layers[1]
    last_layer = layers[end]

    n_in = get_refractive_index(first_layer, λ)
    ε_0in = dielectric_constant(n_in)
    ξ = √(ε_0in) * sin(θ)

    Λ_1324 = @SMatrix [1 0 0 0;
                       0 0 1 0;
                       0 1 0 0;
                       0 0 0 1]
   
    D_0, P_0, γ_0, q_0 = layer_matrices(first_layer, λ, ξ, μ)
    D_f, P_f, γ_f, q_f = layer_matrices(last_layer, λ, ξ, μ)

    
    Γ = I
    Ds = [D_0]
    Ps = Function[P_0]
    γs = [γ_0]
    for layer in layers[2:end - 1]
        D_i, P_i, γ_i, q_i = layer_matrices(layer, λ, ξ, μ)
        T_i = D_i * P_i(layer.thickness) * inv(D_i)
        Γ *= T_i
        push!(Ds, D_i)
        push!(Ps, P_i)
        push!(γs, γ_i)
    end

    push!(Ds, D_f)
    push!(Ps, P_f)
    push!(γs, γ_f)

    Γ = inv(Λ_1324) * inv(D_0) * Γ * D_f * Λ_1324
    r, R, t, T = calculate_tr(Γ)
    S = poynting(ξ, q_0, q_f, γ_0, γ_f, t, r)

    return Γ, S, Ds, Ps, γs
end

"""
    calculate_tr(Γ)

Calculate reflectance and transmittance for the total structure.
This takes the matrix Γ*, but for brevity we call it Γ in this function.

This follows the formalism in:

Yeh, Electromagnetic propagation in birefringent layered media, 1979,
DOI: 10.1364/JOSA.69.000742
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

    r = [rpp, rps, rss, rsp]
    R = [Rpp, Rss, Rsp, Rps]

    t = [tpp, tps, tsp, tss]
    T = [Tpp, Tss, Tsp, Tps]

    return r, R, t, T
end

"""
    calculate_tr(S::Poynting)

Calculate transmittance from the Poynting vector struct,
which contains incident and transmitted energy for both
p-polarized and s-polarized waves.
"""
function calculate_tr(S::Poynting)

    Tpp = S.out_p[3] / S.in_p[3]
    Tss = S.out_s[3] / S.in_s[3]
    
    Rpp = -S.refl_p[3] / S.in_p[3]
    Rss = -S.refl_s[3] / S.in_s[3]

    return Tpp, Tss, Rpp, Rss
end

"""
    calculate_tr(layers, θ=0.0)

Calculate the transmittance and reflectance spectrum
of the structure at a single incidence angle θ.
Accurate transmittance must be calculated via the Poynting
vector. Reflectance is calculated directly from the transfer matrix elements.
"""
function calculate_tr(λ, layers, θ=0.0, μ=1.0+0.0im)
    
    Γ, S, Ds, Ps, γs = propagate(λ, layers, θ, μ)
    r, R, t, T = calculate_tr(Γ)
    Tpp, Tss, Rpp_, Rss_ = calculate_tr(S)
    Rpp = R[1]
    Rss = R[2]
    return Tpp, Tss, Rpp, Rss
end

function angle_resolved(λs, θs, layers)
    Tpp = Matrix{Float64}(undef, length(θs), length(λs))
    Tss = Matrix{Float64}(undef, length(θs), length(λs))
    Rpp = Matrix{Float64}(undef, length(θs), length(λs))
    Rss = Matrix{Float64}(undef, length(θs), length(λs))
    
    Threads.@threads for (i, θ) in collect(enumerate(θs))
        for (j, λ) in collect(enumerate(λs))
            Tpp_, Tss_, Rpp_, Rss_ = calculate_tr(λ, layers, deg2rad(θ))
            Tpp[i, j] = Tpp_
            Tss[i, j] = Tss_
            Rpp[i, j] = Rpp_
            Rss[i, j] = Rss_
        end
    end
    return AngleResolvedResult(Rpp, Rss, Tpp, Tss)
end

function tune_thickness(λs, ts, layers, t_index, θ=0.0)
    Tpp = Matrix{Float64}(undef, length(ts), length(λs))
    Tss = Matrix{Float64}(undef, length(ts), length(λs))
    Rpp = Matrix{Float64}(undef, length(ts), length(λs))
    Rss = Matrix{Float64}(undef, length(ts), length(λs))
    
    Threads.@threads for (i, t) in collect(enumerate(ts))
        changing_layer = Layer(layers[t_index].dispersion, t)
        new_layers = [layers[1:t_index-1]; changing_layer; layers[t_index+1:end]]
        for (j, λ) in collect(enumerate(λs))
            Tpp_, Tss_, Rpp_, Rss_ = calculate_tr(λ, new_layers, θ)
            Tpp[i, j] = Tpp_
            Tss[i, j] = Tss_
            Rpp[i, j] = Rpp_
            Rss[i, j] = Rss_
        end
    end
    return AngleResolvedResult(Rpp, Rss, Tpp, Tss)
end

"""
    electric_field(layers, λ, θ; dz)

Calculate the electric field profile for the entire structure
as a function of z for a given incidence angle θ.
"""
function electric_field(λ, layers, θ=0.0; dz=0.001)

    Γ, S, Ds, Ps, γs = propagate(λ, layers, θ)
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
    Eminus_p[end, :] = inv(P_f(last_layer.thickness)) * Eplus_p[end, :]
    Eminus_s[end, :] = inv(P_f(last_layer.thickness)) * Eplus_s[end, :]

    D_i = Ds[end]

    for l in reverse(eachindex(layers))
        if l >= 2
            layer = layers[l - 1]
            D_prev = Ds[l - 1]
            P_prev = Ps[l - 1]
            L_i = inv(Ds[l - 1]) * D_i

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

    field_p = []
    field_s = []

    field = zeros(ComplexF64, 6, length(zs))

    i = 1
    currentlayer = layers[i]
    for (j, z) in enumerate(zs)
        if z > interface_positions[i]
            i += 1
            currentlayer = layers[i]
        end

        P_i = Ps[i]
        field_p = P_i(-(z - interface_positions[i])) * Eminus_p[i, :]
        field_s = P_i(-(z - interface_positions[i])) * Eminus_s[i, :]

        @views field[1:3, j] = field_p[1] * γs[i][1, :] + field_p[2] * γs[i][2, :] + field_p[3] * γs[i][3, :] + field_p[4] * γs[i][4, :]
        @views field[4:6, j] = field_s[1] * γs[i][1, :] + field_s[2] * γs[i][2, :] + field_s[3] * γs[i][3, :] + field_s[4] * γs[i][4, :]
    end

    return ElectricField(zs, field[1:3, :], field[4:6, :], interface_positions[1:end - 1])
end
