
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
    euler_rotation_matrix(φ, θ, ψ)

Return the 3×3 rotation matrix for ZYZ Euler angles (in radians).

This transforms vectors from the crystal frame to the lab frame:
`v_lab = R * v_crystal`

The rotation is performed as: R = Rz(φ) * Ry(θ) * Rz(ψ)

# Convention
- φ (phi): First rotation about z-axis (0 to 2π)
- θ (theta): Rotation about new y-axis (0 to π) - the tilt angle
- ψ (psi): Second rotation about new z-axis (0 to 2π)

# Common cases
- Optic axis along z: `(0, 0, 0)` - no rotation needed
- Optic axis in xz-plane at angle θ from z: `(0, θ, 0)`
- Quarter-wave plate at 45°: `(π/4, π/2, 0)` for optic axis in xy-plane
"""
function euler_rotation_matrix(φ::Real, θ::Real, ψ::Real)
    cφ, sφ = cos(φ), sin(φ)
    cθ, sθ = cos(θ), sin(θ)
    cψ, sψ = cos(ψ), sin(ψ)

    # ZYZ convention: R = Rz(φ) * Ry(θ) * Rz(ψ)
    @SMatrix [
        cφ*cθ*cψ - sφ*sψ   -cφ*cθ*sψ - sφ*cψ   cφ*sθ;
        sφ*cθ*cψ + cφ*sψ   -sφ*cθ*sψ + cφ*cψ   sφ*sθ;
        -sθ*cψ              sθ*sψ               cθ
    ]
end


"""
    rotate_dielectric_tensor(ε_diag, R)

Rotate a diagonal dielectric tensor from crystal frame to lab frame.

Given a diagonal tensor ε in the crystal's principal axis frame and a rotation
matrix R, returns the rotated tensor: `ε_lab = R * ε * R'`

# Arguments
- `ε_diag`: Diagonal dielectric tensor in crystal frame
- `R`: 3×3 rotation matrix from `euler_rotation_matrix`

# Returns
Full 3×3 SMatrix (may have off-diagonal elements after rotation)
"""
function rotate_dielectric_tensor(ε_diag::Diagonal, R::SMatrix{3,3})
    # ε_lab = R * ε_crystal * R'
    # For diagonal ε, this is equivalent to: R * Diagonal(ε) * R'
    ε_full = SMatrix{3,3,ComplexF64}(ε_diag)
    return R * ε_full * R'
end


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

Supports both isotropic layers (single refractive index) and anisotropic layers
(different refractive indices along principal axes), with optional crystal rotation
via Euler angles.
"""
function layer_matrices(layer, λ, ξ, μ_i)

    ω = 2π * c_0 / λ
    nx, ny, nz = get_refractive_indices(layer, λ)
    εx = dielectric_constant(nx)
    εy = dielectric_constant(ny)
    εz = dielectric_constant(nz)
    ε_diag = dielectric_tensor(εx, εy, εz)

    # Apply rotation if layer has non-zero Euler angles
    φ, θ, ψ = get_euler_angles(layer)
    if φ != 0.0 || θ != 0.0 || ψ != 0.0
        R = euler_rotation_matrix(φ, θ, ψ)
        ε = rotate_dielectric_tensor(ε_diag, R)
    else
        ε = ε_diag
    end

    μ = permeability_tensor(μ_i, μ_i, μ_i)

    M = construct_M(ε, μ)
    a = construct_a(ξ, M)
    Δ = construct_Δ(ξ, M, a)
    q, S = calculate_q(Δ, a)
    q = ComplexF64.(q)
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

# Specialized method for isotropic materials without magnetoelectric coupling (ρ₁ = ρ₂ = 0)
# This is the common case and avoids heap allocations from zeros() and hvcat
function construct_M(ε::Diagonal{ComplexF64,SVector{3,ComplexF64}},
                     μ::Diagonal{ComplexF64,SVector{3,ComplexF64}})
    z = zero(ComplexF64)
    return @SMatrix [
        ε[1,1] z z z z z;
        z ε[2,2] z z z z;
        z z ε[3,3] z z z;
        z z z μ[1,1] z z;
        z z z z μ[2,2] z;
        z z z z z μ[3,3]
    ]
end

# Specialized method for rotated anisotropic materials (full 3×3 dielectric tensor)
function construct_M(ε::SMatrix{3,3,ComplexF64},
                     μ::Diagonal{ComplexF64,SVector{3,ComplexF64}})
    z = zero(ComplexF64)
    return @SMatrix [
        ε[1,1] ε[1,2] ε[1,3] z z z;
        ε[2,1] ε[2,2] ε[2,3] z z z;
        ε[3,1] ε[3,2] ε[3,3] z z z;
        z z z μ[1,1] z z;
        z z z z μ[2,2] z;
        z z z z z μ[3,3]
    ]
end


"""
    construct_a(ξ, M)

Construct the elements of the intermediate 6x6 matrix ``a`` in terms of the
elements of matrix ``M`` (the 6x6 matrix holding the material dielectric and permeability tensors)
and propagation vector ξ. This is implemented as described in 

Berreman, 1972, https://doi.org/10.1364/JOSA.62.000502
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

Berreman, 1972, https://doi.org/10.1364/JOSA.62.000502
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

    return @SMatrix [
        Δ11 Δ12 Δ13 Δ14;
        Δ21 Δ22 Δ23 Δ24;
        Δ31 Δ32 Δ33 Δ34;
        Δ41 Δ42 Δ43 Δ44
    ]
end


"""
    calculate_q(Δ, a)

The four eigenvalues of ``q`` may be obtained from the 
4x4 matrix Δ and then eigenvectors may be found for each eigenvalue.
Here the eigenvalues must be sorted appropriately to avoid 
potentially discontinuous solutions. This extends from the work in

Li et al, 1988, https://doi.org/10.1364/AO.27.001334
"""
function calculate_q(Δ, a)

    q_unsorted, Ψ_unsorted = eigen(Δ)
    transmitted_mode = MVector(0, 0)
    reflected_mode = MVector(0, 0)

    kt = 1
    kr = 1
    if isreal(q_unsorted)
        for m in 1:4
            if real(q_unsorted[m]) >= 0.0
                kt <= 2 && (transmitted_mode[kt] = m)
                kt += 1
            else
                kr <= 2 && (reflected_mode[kr] = m)
                kr += 1
            end
        end
    else
        for m in 1:4
            if imag(q_unsorted[m]) >= 0.0
                kt <= 2 && (transmitted_mode[kt] = m)
                kt += 1
            else
                kr <= 2 && (reflected_mode[kr] = m)
                kr += 1
            end
        end
    end

    # Verify invariant: exactly 2 transmitted and 2 reflected modes
    if kt != 3 || kr != 3
        throw(ArgumentError(
            "Mode sorting failed: expected 2 transmitted and 2 reflected modes, " *
            "got $(kt-1) transmitted and $(kr-1) reflected. Eigenvalues: $q_unsorted"
        ))
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

This is based on the work in
Xu, et al., 2000, https://doi.org/10.1103/PhysRevB.61.1740
"""
function calculate_γ(ξ, q, ε, μ)

    γ = @MMatrix zeros(ComplexF64, 4, 3)
    γ[1,1] = 1
    γ[2,2] = 1
    γ[4,2] = 1
    γ[3,1] = -1

    # Common denominator term - can become singular at specific angles
    denom_33 = μ * ε[3,3] - ξ^2
    singular_33 = abs(denom_33) < eps(Float64)

    if isapprox(q[1], q[2])
        γ[1,2] = 0
        γ[2,1] = 0
    else
        γ[1,2] = (μ * ε[2,3] * (μ * ε[3,1] + ξ * q[1]) - μ * ε[2,1] * denom_33) / (denom_33 * (μ * ε[2,2] - ξ^2 - q[1]^2) - μ^2 * ε[2,3] * ε[3,2])
        γ[2,1] = (μ * ε[3,2] * (μ * ε[1,3] + ξ * q[2]) - μ * ε[1,2] * denom_33) / (denom_33 * (μ * ε[1,1] - q[2]^2) - (μ * ε[1,3] + ξ * q[2]) * (μ * ε[3,1] + ξ * q[2]))
    end

    # γ[i,3] components use denom_33 directly - set to zero if singular
    if singular_33
        γ[1,3] = 0
        γ[2,3] = 0
    else
        γ[1,3] = (-μ * ε[3,1] - ξ * q[1] - μ * ε[3,2] * γ[1,2]) / denom_33
        γ[2,3] = (-(μ * ε[3,1] + ξ * q[2]) * γ[2,1] - μ * ε[3,2]) / denom_33
    end

    if isapprox(q[3], q[4])
        γ[3,2] = 0.0
        γ[4,1] = 0.0
    else
        γ[3,2] = (μ * ε[2,1] * denom_33 - μ * ε[2,3] * (μ * ε[3,1] + ξ * q[3])) / (denom_33 * (μ * ε[2,2] - ξ^2 - q[3]^2) - μ^2 * ε[2,3] * ε[3,2])
        γ[4,1] = (μ * ε[3,2] * (μ * ε[1,3] + ξ * q[4]) - μ * ε[1,2] * denom_33) / (denom_33 * (μ * ε[1,1] - q[4]^2) - (μ * ε[1,3] + ξ * q[4]) * (μ * ε[3,1] + ξ * q[4]))
    end

    # γ[i,3] components for backward modes - set to zero if singular
    if singular_33
        γ[3,3] = 0
        γ[4,3] = 0
    else
        γ[3,3] = (μ * ε[3,1] + ξ * q[3] + μ * ε[3,2] * γ[3,2]) / denom_33
        γ[4,3] = (-(μ * ε[3,1] + ξ * q[4]) * γ[4,1] - μ * ε[3,2] ) / denom_33
    end

    # normalize γ (use SVector to avoid slice allocations)
    # Skip normalization if the row is effectively zero to avoid NaN from 0/0
    for i in 1:4
        v = SVector(γ[i,1], γ[i,2], γ[i,3])
        Z = √(v ⋅ v')
        if abs(Z) > eps(Float64)
            γ[i,1] /= Z
            γ[i,2] /= Z
            γ[i,3] /= Z
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

Xu, et al., 2000, https://doi.org/10.1103/PhysRevB.61.1740
"""
function dynamical_matrix(ξ, q, γ, μ)
    A_3 = (γ[:, 1] .* q .- ξ * γ[:, 3]) ./ μ
    A_4 = γ[:, 2] .* q ./ μ

    return @SMatrix [
        γ[1, 1] γ[2, 1] γ[3, 1] γ[4, 1];
        γ[1, 2] γ[2, 2] γ[3, 2] γ[4, 2];
        A_3[1] A_3[2] A_3[3] A_3[4];
        A_4[1] A_4[2] A_4[3] A_4[4]
    ]
end


struct PropagationMatrix{Tω,Tq}
    ω::Tω
    q::Tq
end

"""
    propagation_matrix(ω, q)

Returns a callable object that propagates the electromagnetic
field a distance z through a material for a frequency ω
and wavevector `q`.

The propagation matrix is diagonal with elements:
```
P(z) = diag(exp(-iωq₁z/c), exp(-iωq₂z/c), exp(-iωq₃z/c), exp(-iωq₄z/c))
```

This uses the **exp(-iωt)** time convention, consistent with Berreman (1972)
and Passler & Paarmann (2017). Note that Yeh uses exp(+iωt), which would
flip the sign in the exponent.
"""
propagation_matrix(ω, q) = PropagationMatrix(ω, q)

(P::PropagationMatrix)(z) = Diagonal(exp.(-im * P.ω * P.q * z / c_0))
