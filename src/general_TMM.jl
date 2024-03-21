"""
    construct_a(ξ, M)

Construct the elements of the intermediate 6x6 matrix ``a`` in terms of the
elements of matrix ``M`` (the 6x6 matrix holding the material dielectric and permeability tensors)
and propagation vector ξ. This is implemented as described in 

Berreman, 1972. https://doi.org/10.1364/JOSA.62.000502
"""
function construct_a(ξ, M)

    d = M[3,3] * M[6,6] - M[3,6] * M[6,3]

    a31 = (M[6,1] * M[3,6] - M[3,1] * M[6,6]) / d
    a32 =((M[6,2] - ξ) * M[3,6] - M[3,2] * M[6,6]) / d
    a34 = (M[6,4] * M[3,6] -  M[3,4] * M[6,6]) / d
    a35 = (M[6,5] * M[3,6] - (M[3,5] + ξ) * M[6,6]) / d

    a61 = (M[6,3] * M[3,1] - M[3,3] * M[6,1]) / d
    a62 = (M[6,3] * M[3,2] - M[3,3] * (M[6,2] - ξ)) / d
    a64 = (M[6,3] * M[3,4] - M[3,3] * M[6,4]) / d
    a65 = (M[6,3] * (M[3,5] + ξ) - M[3,3] * M[6,5]) / d
    
    return @SMatrix [
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        a31 a32 0 a34 a35 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        a61 a62 0 a64 a65 0
    ]
end

"""
    construct_Δ(ξ, M, a)

Construct the reordered matrix Δ in terms of the elements of
the two matrices, M and a, and the in-plane reduced wavevector ξ = ``k_x / k_0``.
The matrix Δ is involved in the relation

```math
    \\frac{\\delta}{\\delta z}\\Psi = \\frac{i \\omega}{c}\\Delta \\Psi
```

and Δ is the reordered S matrix in Berreman's formulation.

Berreman, 1972. https://doi.org/10.1364/JOSA.62.000502
"""
function construct_Δ(ξ, M, a)

    Δ = (i, j) -> M[i,j] + M[i,3] * a[3,j] + (i == 2 || i == 4 ? M[i,6] - ξ : M[i,6]) * a[6,j]

    return @SMatrix [
        Δ(5,1) Δ(5,5) Δ(5,2) Δ(5,4);
        Δ(1,1) Δ(1,5) Δ(1,2) Δ(1,4);
        Δ(4,1) Δ(4,5) Δ(4,2) Δ(4,4);
        Δ(2,1) Δ(2,5) Δ(2,2) Δ(2,4)
    ]
end

"""
    calculate_q(Δ, a)

The four eigenvalues of ``q`` may be obtained from the 
4x4 matrix Δ and then eigenvectors may be found for each eigenvalue.
Here the eigenvalues must be sorted appropriately to avoid 
potentially discontinuous solutions. This extends from the work in

Li et al, 1988. https://doi.org/10.1364/AO.27.001334
"""
function calculate_q(Δ, a)

    q_unsorted, Ψ_unsorted = eigen(Δ)

    transmitted_mode = @MVector zeros(Int, 2)
    reflected_mode = @MVector zeros(Int, 2)

    kt = 1
    kr = 1

    if isreal(q_unsorted)
        for m in 1:4
            if real(q_unsorted[m]) >= 0.0
                transmitted_mode[kt] = m
                kt += 1
            else
                reflected_mode[kr] = m
                kr += 1
            end
        end
    else
        for m in 1:4
            if imag(q_unsorted[m]) >= 0.0
                transmitted_mode[kt] = m
                kt += 1
            else
                reflected_mode[kr] = m
                kr += 1
            end
        end
    end

    S_unsorted = poynting(Ψ_unsorted, a)
    t, r = evaluate_birefringence(Ψ_unsorted, S_unsorted, transmitted_mode, reflected_mode)

    q = SVector(q_unsorted[t[1]], q_unsorted[t[2]], q_unsorted[r[1]], q_unsorted[r[2]])
    S = SVector(S_unsorted[t[1]], S_unsorted[t[2]], S_unsorted[r[1]], S_unsorted[r[2]])

    return q, S
end

"""
    calculate_γ(ξ, q, ε, μ)

The 4 x 3 matrix γ contains vector components that belong 
to the electric field calculated such that singularities are identified and removed.

`q[1]` and `q[2]` are forward-traveling modes and
`q[3]` and `q[4]` are backward-traveling modes.

This is based on the work in:

Xu et al. Optical degeneracies in anisotropic layered media:
Treatment of singularities in a 4x4 matrix formalism, 2000.
DOI: 10.1103/PhysRevB.61.1740
"""
function calculate_γ(ξ, q, ε, μ)

    γ = @MMatrix zeros(ComplexF64, 4, 3)

    γ[1,1] = 1
    γ[2,2] = 1
    γ[4,2] = 1
    γ[3,1] = -1

    if isapprox(q[1], q[2])
        γ[1,2] = 0
        γ[2,1] = 0
    else
        γ[1,2] = (μ * ε[2,3] * (μ * ε[3,1] + ξ * q[1]) - μ * ε[2,1] * (μ * ε[3,3] - ξ^2)) / ((μ * ε[3,3] - ξ^2) * (μ * ε[2,2] - ξ^2 - q[1]^2) - μ^2 * ε[2,3] * ε[3,2])
        γ[2,1] = (μ * ε[3,2] * (μ * ε[1,3] + ξ * q[2]) - μ * ε[1,2] * (μ * ε[3,3] - ξ^2)) / ((μ * ε[3,3] - ξ^2) * (μ * ε[1,1] - q[2]^2) - (μ * ε[1,3] + ξ * q[2]) * (μ * ε[3,1] + ξ * q[2]))
    end

    γ[1,3] = (-μ * ε[3,1] - ξ * q[1] - μ * ε[3,2] * γ[1,2]) / (μ * ε[3,3] - ξ^2)
    γ[2,3] = (-(μ * ε[3,1] + ξ * q[2]) * γ[2,1] - μ * ε[3,2]) / (μ * ε[3,3] - ξ^2)

    if isapprox(q[3], q[4])
        γ[3,2] = 0.0
        γ[4,1] = 0.0
    else
        γ[3,2] = (μ * ε[2,1] * (μ * ε[3,3] - ξ^2) - μ * ε[2,3] * (μ * ε[3,1] + ξ * q[3])) / ((μ * ε[3,3] - ξ^2) * (μ * ε[2,2] - ξ^2 - q[3]^2) - μ^2 * ε[2,3] * ε[3,2])
        γ[4,1] = (μ * ε[3,2] * (μ * ε[1,3] + ξ * q[4]) - μ * ε[1,2] * (μ * ε[3,3] - ξ^2)) / ((μ * ε[3,3] - ξ^2) * (μ * ε[1,1] - q[4]^2) - (μ * ε[1,3] + ξ * q[4]) * (μ * ε[3,1] + ξ * q[4]))
    end

    γ[3,3] = (μ * ε[3,1] + ξ * q[3] + μ * ε[3,2] * γ[3,2]) / (μ * ε[3,3] - ξ^2)
    γ[4,3] = (-(μ * ε[3,1] + ξ * q[4]) * γ[4,1] - μ * ε[3,2] ) / (μ * ε[3,3] - ξ^2)

    # normalize γ
    for i in 1:4
        Z = √(γ[i, :] ⋅ γ[i, :]')
        for j in 1:3
            γ[i, j] /= Z
        end
    end

    return SMatrix(γ)
end

"""
    dynamical_matrix(ξ, q, γ, μ)

The dynamical matrix relating two layers at the interface
where matrix ``A_i`` for layer ``i`` relates the field ``E_i`` to
the field in the previous layer ``i - 1`` via

```math
    A_{i-1}E_{i-1} = A_{i}E_{i}
```

Xu et al., 2000,
DOI: 10.1103/PhysRevB.61.1740
"""
dynamical_matrix(ξ, q, γ, μ) = @SMatrix{4, 4}(
        γ[:, 1],
        γ[:, 2],
        (γ[:, 1] .* q .- ξ * γ[:, 3]) ./ μ,
        γ[:, 2] .* q ./ μ
    )

"""
    propagation_matrix(ω, q)

Returns a function that propagates the electromagnetic
field a distance z through a material for a frequency ω
and wavevector ``q``.
"""
function propagation_matrix(ω, q)
    return z -> Diagonal(exp.(-im * ω * q * z / c_0))
end

"""
    transfermatrix(ω, ξ, q, γ, μ, d)

Calcuate the transfer matrix from the dynamical (interface) matrix
and propagation matrix for one layer.
"""
function transfermatrix(ω, ξ, q, γ, μ, d)

    A_i = dynamical_matrix(ξ, q, γ, μ)
    P_i = propagation_matrix(ω, q)
    T_i = A_i * P_i(d) * inv(A_i)

    return A_i, P_i, T_i
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
    
    # Note: isapprox(NaN, NaN) = false, which is important
    # in the case that both Sx and Sy are zero.
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
function abs_ratio(a, b)
    return abs2(a) / (abs2(a) + abs2(b))
end

"""
    electric_field(s::Structure, λ, θ; numpoints)

Calculate the electric field profile for the entire structure
as a function of z for a given incidence angle θ.
"""
function electric_field(s::Structure, λ, θ = 0.0; numpoints = 1000)

    s = Structure(s.layers, [λ], [θ])
    res = calculate_Γ_S(s, θ)
    rs, Rs, ts, Ts = tr_from_Γ(res.tm)

    ω = 2π * c_0 / s.λ[1]
    μ = 1.0 + 0.0im

    t = ts[1]
    ξ = res.ξ[1]

    superstrate = s.layers[1]
    substrate = s.layers[end]
    A_f, P_f, T_f, γ_f, q_f = layer_params(ω, ξ, substrate.n_r[1] + substrate.n_i[1] * im, μ, substrate.thickness)

    Eplus_p = zeros(ComplexF64, length(s.layers), 4)
    Eminus_p = zeros(ComplexF64, length(s.layers), 4)

    Eplus_s = zeros(ComplexF64, length(s.layers), 4)
    Eminus_s = zeros(ComplexF64, length(s.layers), 4)

    Eplus_p[end, :] = [t[1], t[2], 0, 0]
    Eplus_s[end, :] = [t[3], t[4], 0, 0]
    
    Eminus_p[end, :] = inv(P_f(substrate.thickness)) * Eplus_p[end, :]
    Eminus_s[end, :] = inv(P_f(substrate.thickness)) * Eplus_s[end, :]

    propagation_funcs = Function[P_f]
    γs = [γ_f]
    A_i = A_f

    for l in reverse(eachindex(s.layers))
    
        if l >= 2

            layer = s.layers[l - 1]
            
            n = layer.n_r[1] + layer.n_i[1] * im
            A_prev, P_prev, T_prev, γ_prev, q_prev = layer_params(ω, ξ, n, μ, layer.thickness)

            push!(propagation_funcs, P_prev)
            push!(γs, γ_prev)

            L_i = inv(A_prev) * A_i

            Eminus_p[l - 1, :] = L_i * Eplus_p[l, :]
            Eminus_s[l - 1, :] = L_i * Eplus_s[l, :]

            Eplus_p[l - 1, :] = P_prev(layer.thickness) * Eminus_p[l - 1, :]
            Eplus_s[l - 1, :] = P_prev(layer.thickness) * Eminus_s[l - 1, :]
            
            A_i = A_prev
        end
    end

    reverse!(propagation_funcs)
    reverse!(γs)

    interface_positions, total_thickness = find_layerbounds(s)
    interface_positions .-= superstrate.thickness

    zs = range(-superstrate.thickness, interface_positions[end], length = numpoints)

    field_p = []
    field_s = []

    field = zeros(ComplexF64, 6, length(zs))

    i = 1
    currentlayer = s.layers[i]
    for (j, z) in enumerate(zs)
        
        if z > interface_positions[i]
            i += 1
            currentlayer = s.layers[i]
        end

        P_i = propagation_funcs[i]
        field_p = P_i(-(z - interface_positions[i])) * Eminus_p[i, :]
        field_s = P_i(-(z - interface_positions[i])) * Eminus_s[i, :]

        field[1:3, j] = field_p[1] * γs[i][1, :] + field_p[2] * γs[i][2, :] + field_p[3] * γs[i][3, :] + field_p[4] * γs[i][4, :]
        field[4:6, j] = field_s[1] * γs[i][1, :] + field_s[2] * γs[i][2, :] + field_s[3] * γs[i][3, :] + field_s[4] * γs[i][4, :]
     
    end

    return ElectricField(zs, field[1:3, :], field[4:6, :], interface_positions[1:end - 1])
end


"""
    layer_params(ω, ξ, n, μ, d)

Calculate all parameters for a single layer, particularly
the propagation matrix and dynamical matrix so that
the overall transfer matrix can be calculated.
"""
function layer_params(ω, ξ, n, μ, d)

    ε_i = dielectric_constant(n)
    ε = dielectric_tensor(ε_i, ε_i, ε_i)

    M = @MMatrix zeros(ComplexF64, 6, 6)
    M[1:3, 1:3] = ε
    M[4:6, 4:6] = Diagonal([μ, μ, μ])
    M = SMatrix(M)
    
    a = construct_a(ξ, M)
    Δ = construct_Δ(ξ, M, a)
    q, S = calculate_q(Δ, a)
    γ = calculate_γ(ξ, q, ε, μ)
    A, P, T = transfermatrix(ω, ξ, q, γ, μ, d)

    return A, P, T, γ, q
end

"""
    calculate_Γ_S(s::Structure, θ=0.0)

Calculate the total structure transfer matrix Γ
and the Poynting vector S for all
input frequencies ω, returning the TransferMatrixResult
type, which contains all transfer matrices, Poynting vectors,
and ξ.
"""
function calculate_Γ_S(s::Structure, θ=0.0)

    superstrate = s.layers[1]
    substrate = s.layers[end]

    μ = 1.0 + 0.0im
    ε_0in = dielectric_constant(superstrate)
    ξs = @. √(ε_0in) * sin(θ)

    Λ_1324 = @SMatrix [1 0 0 0;
                      0 0 1 0;
                      0 1 0 0;
                      0 0 0 1]

    Γs = AbstractVector{SMatrix{4, 4, ComplexF64}}([])
    Ss = AbstractVector{Poynting}([])

    ωs = 2π * c_0 ./ s.λ

    for (i, ω) in enumerate(ωs)
        
        A_0, P_0, T_0, γ_0, q_0 = layer_params(ω, ξs[i], superstrate.n_r[i] + superstrate.n_i[i] * im, μ, superstrate.thickness)
        A_f, P_f, T_f, γ_f, q_f = layer_params(ω, ξs[i], substrate.n_r[i] + substrate.n_i[i] * im, μ, substrate.thickness)
        
        Γ = I

        for layer in s.layers[2:end - 1]

            n = layer.n_r[i] + layer.n_i[i] * im
            A_i, P_i, T_i, γ_i, q_i = layer_params(ω, ξs[i], n, μ, layer.thickness)
            Γ *= T_i
        end

        Γ = inv(Λ_1324) * inv(A_0) * Γ * A_f * Λ_1324
        r, R, t, T = tr_from_Γ(Γ)
        S = poynting(ξs[i], q_0, q_f, γ_0, γ_f, t, r)
        push!(Γs, Γ)
        push!(Ss, S)
    end

    return TransferMatrixResult(Γs, Ss, ξs)
end

"""
    tr_from_Γ(Γ::Matrix)

Calculate reflectance and transmittance for the total structure.
This takes the matrix Γ*, but for brevity we call it Γ in this function.

This follows the formalism in:

Yeh, Electromagnetic propagation in birefringent layered media, 1979,
DOI: 10.1364/JOSA.69.000742
"""
function tr_from_Γ(Γ::SMatrix)

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
    tr_from_Γ(Γs::Vector)

Calculate reflectance and transmittance for the total structure for
each `Γ` in the Vector of `Γs`.
"""
function tr_from_Γ(Γs::Vector)

    rs = []
    Rs = []
    ts = []
    Ts = []

    for Γ in Γs
        r, R, t, T = tr_from_Γ(Γ)
        push!(rs, r)
        push!(Rs, R)
        push!(ts, t)
        push!(Ts, T)
    end

    return rs, Rs, ts, Ts
end

"""
    tr_from_poynting(S::Poynting)

Calculate transmittance from the Poynting vector struct,
which contains incident and transmitted energy for both
p-polarized and s-polarized waves.
"""
function tr_from_poynting(S::Poynting)

    Tpp = S.out_p[3] / S.in_p[3]
    Tss = S.out_s[3] / S.in_s[3]

    Rpp = - S.refl_p[3] / S.in_p[3]
    Rss = - S.refl_s[3] / S.in_s[3]

    return Tpp, Tss, Rpp, Rss
end

"""
    tr_from_poynting(Ss::Vector{Poynting})

Calculate transmission and reflection spectra from a Vector
of `Poynting` types.
"""
function tr_from_poynting(Ss::Vector{Poynting})

    Tpps = Float64[]
    Tsss = Float64[]
    Rpps = Float64[]
    Rsss = Float64[]

    for S in Ss
        Tpp, Tss, Rpp, Rss = tr_from_poynting(S)
        push!(Tpps, Tpp)
        push!(Tsss, Tss)
        push!(Rpps, Rpp)
        push!(Rsss, Rss)
    end

    return Tpps, Tsss, Rpps, Rsss
end

"""
    angle_resolved(s::Structure)

Iterate through each angle provided in the structure
to find the reflectance and transmittance spectra from
the calculated transfer matrices and Poynting vectors.
"""
function angle_resolved(s::Structure)

    Rpp_spectrum = Matrix{Real}(undef, length(s.θ), length(s.λ))
    Rss_spectrum = Matrix{Real}(undef, length(s.θ), length(s.λ))
    Tpp_spectrum = Matrix{Real}(undef, length(s.θ), length(s.λ))
    Tss_spectrum = Matrix{Real}(undef, length(s.θ), length(s.λ))
    Γs = []
    ξ = Matrix{ComplexF64}(undef, length(s.θ), length(s.λ))

    for (i, θ) in enumerate(s.θ)

        result = calculate_Γ_S(s, θ)

        rs, Rs, ts, Ts = tr_from_Γ(result.tm)

        ξ[i, :] = result.ξ

        Tpp, Tss, Rpp, Rss = tr_from_poynting(result.poynting)
        Rpp_spectrum[i, :] = [R[1] for R in Rs]
        Rss_spectrum[i, :] = [R[2] for R in Rs]
        # Rpp_spectrum[i, :] = Rpp
        # Rss_spectrum[i, :] = Rss
        Tpp_spectrum[i, :] = Tpp
        Tss_spectrum[i, :] = Tss
    end

    return AngleResolvedResult(Rpp_spectrum, Rss_spectrum, Tpp_spectrum, Tss_spectrum, Γs, ξ)
end

"""
    calculate_tr(s::Structure, θ=0.0)

Calculate the transmittance and reflectance spectrum
of the structure at a single incidence angle θ.
Accurate transmittance must be calculated via the Poynting
vector. Reflectance is calculated directly from the transfer matrix elements.
"""
function calculate_tr(s::Structure, θ=0.0)
    
    result = calculate_Γ_S(s, θ)
    rs, Rs, ts, Ts = tr_from_Γ(result.tm)
    ξ = result.ξ
    Tpp, Tss, Rpp, Rss = tr_from_poynting(result.poynting)

    Rpp = [R[1] for R in Rs]
    Rss = [R[2] for R in Rs]

    return Tpp, Tss, Rpp, Rss
end

"""
    initialize(structure::Structure, λs)

Initializing a `Structure` interpolates the wavelength-dependent
refractive index data using the given `λs` Vector for all `Layer`s
in the `Structure`, returning a new structure with the interpolated data.
"""
function initialize(structure::Structure, λs)
    layers = Layer[]
    for layer in structure.layers
        new_layer = interp_data(layer, λs)
        push!(layers, new_layer)
    end
    return Structure(layers, collect(λs), structure.θ)
end

"""
    interp_data(layer::Layer, λs)

Given a new set of wavelengths (the Vector, λs), interpolate the 
complex refractive index values for the input Layer
and return a new Layer with the new λ, n_r, and n_i.
The new wavelengths must not extend beyond the 
domain of the existing wavelengths in the Layer (i.e. no extrapolation).

Here we use LinearInterpolation from the package DataInterpolations.jl
"""
function interp_data(layer::Layer, λs)

    if length(layer.λ) == 1
        n_r = fill(layer.n_r[1], length(λs))
        n_i = fill(layer.n_i[1], length(λs))

        return Layer(layer.material, layer.thickness, λs, n_r, n_i)
    else
        interp_n_r = LinearInterpolation(layer.n, layer.λ)
        interp_n_i = LinearInterpolation(layer.n_i, layer.λ)

        n_r = interp_n_r.(λs)
        n_i = interp_n_i.(λs)

        return Layer(layer.material, layer.thickness, λs, n_r, n_i)
    end
end

"""
    find_layerbounds(s::Structure)

Find the unitful z coordinate for all layer-layer interfaces in the structure,
with the first interface starting at z = 0.
(negative z corresponds to positions inside the first layer.)
"""
function find_layerbounds(s::Structure)

    total_thickness = 0.0
    interface_positions = Float64[]
    
    for layer in s.layers
        push!(interface_positions, total_thickness + layer.thickness)
        total_thickness += layer.thickness
    end

    return interface_positions, total_thickness
end
