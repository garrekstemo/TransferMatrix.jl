
"""
    dielectric_tensor(ε1, ε2, ε3)

Return the diagonal complex dielectric tensor

```math
\\varepsilon = 
\\begin{pmatrix}
    \\varepsilon_1 & 0 & 0 \\\\\
    0 & \\varepsilon_2  & 0 \\\\\
    0 & 0 & \\varepsilon_3
\\end{pmatrix}
```
"""
dielectric_tensor(ε1, ε2, ε3) = Diagonal(SVector{3, ComplexF64}(ε1, ε2, ε3))


"""
    permeability_tensor(μ1, μ2, μ3)

This produces the diagonal permeability tensor, 
which is identical to the way we build the `dielectric_tensor`,
and we include this function simply for completeness.
"""
permeability_tensor(μ1, μ2, μ3) = Diagonal(SVector{3, ComplexF64}(μ1, μ2, μ3))


"""
    layer_matrices(ω, ξ, layer, μ)

Calculate all parameters for a single layer, particularly
the propagation matrix and dynamical matrix so that
the overall transfer matrix can be calculated.
"""
function layer_matrices(layer, λ, ξ, μ_i)

    ω = 2π * c_0 / λ
    n_i = get_refractive_index(layer.material, λ)
    ε_i = dielectric_constant(n_i)
    ε = dielectric_tensor(ε_i, ε_i, ε_i)
    μ = permeability_tensor(μ_i, μ_i, μ_i)

    M = construct_M(ε, μ)
    a = construct_a(ξ, M)
    Δ = construct_Δ(ξ, M, a)
    q, S = calculate_q(Δ, a)
    γ = calculate_γ(ξ, q, ε, μ_i)
    D = dynamical_matrix(ξ, q, γ, μ_i)
    P = propagation_matrix(ω, q)
    return D, P, γ, q
end

"""
    construct_M(ε, μ, ρ1, ρ2)

Construct the 6x6 matrix M from the dielectric and permeability tensors.
"""
function construct_M(ε, μ=Diagonal(ones(3)), ρ1=zeros(3, 3), ρ2=zeros(3, 3))
    return [ε ρ1; ρ2 μ]
end


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
    
    return SMatrix{6, 6}([
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        a31 a32 0 a34 a35 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        a61 a62 0 a64 a65 0
    ])
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

    Δ11 =  M[5,1] + (M[5,3] + ξ) * a[3,1] + M[5,6] * a[6,1]
    Δ12 =  M[5,5] + (M[5,3] + ξ) * a[3,5] + M[5,6] * a[6,5]
    Δ13 =  M[5,2] + (M[5,3] + ξ) * a[3,2] + M[5,6] * a[6,2]
    Δ14 = -M[5,4] - (M[5,3] + ξ) * a[3,4] - M[5,6] * a[6,4]

    Δ21 =  M[1,1] + M[1,3] * a[3,1] + M[1,6] * a[6,1]
    Δ22 =  M[1,5] + M[1,3] * a[3,5] + M[1,6] * a[6,5]
    Δ23 =  M[1,2] + M[1,3] * a[3,2] + M[1,6] * a[6,2]
    Δ24 = -M[1,4] - M[1,3] * a[3,4] - M[1,6] * a[6,4]

    Δ31 = -M[4,1] - M[4,3] * a[3,1] - M[4,6] * a[6,1]
    Δ32 = -M[4,5] - M[4,3] * a[3,5] - M[4,6] * a[6,5]
    Δ33 = -M[4,2] - M[4,3] * a[3,2] - M[4,6] * a[6,2]
    Δ34 =  M[4,4] + M[4,3] * a[3,4] + M[4,6] * a[6,4]

    Δ41 =  M[2,1] + M[2,3] * a[3,1] + (M[2,6] - ξ) * a[6,1]
    Δ42 =  M[2,5] + M[2,3] * a[3,5] + (M[2,6] - ξ) * a[6,5]
    Δ43 =  M[2,2] + M[2,3] * a[3,2] + (M[2,6] - ξ) * a[6,2]
    Δ44 = -M[2,4] - M[2,3] * a[3,4] - (M[2,6] - ξ) * a[6,4]

    return SMatrix{4, 4}([
        Δ11 Δ12 Δ13 Δ14;
        Δ21 Δ22 Δ23 Δ24;
        Δ31 Δ32 Δ33 Δ34;
        Δ41 Δ42 Δ43 Δ44
    ])
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
    transmitted_mode = [0, 0]
    reflected_mode = [0, 0]

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
function dynamical_matrix(ξ, q, γ, μ)
    A_3 = (γ[:, 1] .* q .- ξ * γ[:, 3]) ./ μ
    A_4 = γ[:, 2] .* q ./ μ

    return SMatrix{4, 4}([
        γ[1, 1] γ[2, 1] γ[3, 1] γ[4, 1];
        γ[1, 2] γ[2, 2] γ[3, 2] γ[4, 2];
        A_3[1] A_3[2] A_3[3] A_3[4];
        A_4[1] A_4[2] A_4[3] A_4[4]
    ])
end


"""
    propagation_matrix(ω, q)

Returns a function that propagates the electromagnetic
field a distance z through a material for a frequency ω
and wavevector ``q``.
"""
function propagation_matrix(ω, q)
    return z -> Diagonal(exp.(-im * ω * q * z / c_0))
end