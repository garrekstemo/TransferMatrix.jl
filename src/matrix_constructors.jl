
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
    gyrotropic_tensor(d, od; axis=:z)

Constant gyrotropic (e.g. gyromagnetic Polder) tensor with diagonal `d` and
antisymmetric imaginary off-diagonal `±i·od`, with the gyration vector along
`axis`. For `axis=:z`:

```
[ d      i·od   0
 −i·od   d      0
  0      0      1 ]
```

Hermitian (hence lossless) for real `d, od`. The sign/handedness follows the
package `exp(-iωt)` convention (validated against the `t_ps = -t_sp`
non-reciprocity relation).
"""
function gyrotropic_tensor(d, od; axis::Symbol=:z)
    dC = ComplexF64(d); g = im * ComplexF64(od); o = one(ComplexF64); z = zero(ComplexF64)
    if axis === :z
        return @SMatrix [dC g z; -g dC z; z z o]
    elseif axis === :x
        return @SMatrix [o z z; z dC g; z -g dC]
    elseif axis === :y
        return @SMatrix [dC z g; z o z; -g z dC]
    else
        throw(ArgumentError("axis must be :x, :y, or :z, got :$axis"))
    end
end


"""
    polder_permeability(; f0, fm, linewidth=0.0, axis=:z)

Return a function `f -> μ_tensor` giving the gyromagnetic (Polder) permeability
tensor of a saturated ferrite at frequency `f`:

```math
μ(f) = 1 + \\frac{f_0 f_m}{f_0^2 - f^2}, \\qquad κ(f) = \\frac{f f_m}{f_0^2 - f^2}
```

with the gyration axis `axis`. `f0` is the ferromagnetic-resonance frequency
(`= γμ₀H₀/2π`), `f_m = γμ₀M_s/2π`, and `γ = g·e/2mₑ`. A nonzero `linewidth`
(ΔH, in the same units as `f0`) adds loss via `f0 → f0 − i·linewidth/2`.

Use as a layer permeability with `Layer(...; mu = λ -> polder_permeability(...)(f(λ)))`,
mapping wavelength to frequency as appropriate for your unit system.
"""
function polder_permeability(; f0, fm, linewidth=0.0, axis::Symbol=:z)
    f0c = ComplexF64(f0) - im*ComplexF64(linewidth)/2
    fmC = ComplexF64(fm)
    return function (f)
        denom = f0c^2 - ComplexF64(f)^2
        μ = 1 + f0c*fmC/denom
        κ = ComplexF64(f)*fmC/denom
        return gyrotropic_tensor(μ, κ; axis=axis)
    end
end


# Field-slot reorder swapping slots 2↔3 (transposition (2 3)); its own inverse
# (swap23² = I), so `swap23 \ X` and `swap23 * X` coincide.
const _swap23 = @SMatrix [1.0 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

# Rotation-aware dielectric tensor for a layer at wavelength λ. Shared by the
# eigenmode path (`layer_matrices`) and the matrix-exponential path
# (`layer_transfer_exp`) so ε construction lives in one place.
function _layer_epsilon(layer, λ)
    nx, ny, nz = get_refractive_indices(layer, λ)
    ε_diag = dielectric_tensor(dielectric_constant(nx), dielectric_constant(ny), dielectric_constant(nz))
    φ, θ, ψ = get_euler_angles(layer)
    return (φ != 0.0 || θ != 0.0 || ψ != 0.0) ?
        rotate_dielectric_tensor(ε_diag, euler_rotation_matrix(φ, θ, ψ)) : ε_diag
end

"""
    layer_matrices(ω, k_par, layer, μ)

Calculate all parameters for a single layer, particularly
the propagation matrix and dynamical matrix so that
the overall transfer matrix can be calculated.

Supports both isotropic layers (single refractive index) and anisotropic layers
(different refractive indices along principal axes), with optional crystal rotation
via Euler angles.
"""
function layer_matrices(layer, λ, k_par, μ_i)
    ω = 2π * c_0 / λ
    ε = _layer_epsilon(layer, λ)

    if ismagnetic(layer)
        μmat = get_permeability(layer, λ)                 # SMatrix{3,3,ComplexF64}
        C = construct_constitutive(ε, μmat)  # 6×6 constitutive matrix [ε ρ1; ρ2 μ]
        a = construct_a(k_par, C); Δ = construct_Δ(k_par, C, a)
        q, S = calculate_q(Δ, a); q = ComplexF64.(q)
        E_modes = calculate_E_modes_tensor(k_par, q, ε, μmat)
        D = dynamical_matrix(k_par, q, E_modes, μmat)
    else
        μ = permeability_tensor(μ_i, μ_i, μ_i)
        C = construct_constitutive(ε, μ)
        a = construct_a(k_par, C); Δ = construct_Δ(k_par, C, a)
        q, S = calculate_q(Δ, a); q = ComplexF64.(q)
        E_modes = calculate_E_modes(k_par, q, ε, μ_i)
        D = dynamical_matrix(k_par, q, E_modes, μ_i)
    end
    P = propagation_matrix(ω, q)
    return D, P, E_modes, q
end

"""
    layer_transfer_exp(layer, λ, k_par, ω, μ_i)

Interior-layer 4×4 transfer matrix in the dynamical-matrix field basis
`(Eₓ, Eᵧ, Hᵧ, −Hₓ)`, computed as the **matrix exponential** of the Berreman Δ
matrix rather than by eigenmode decomposition:

```math
T = swap23\\, \\exp\\!\\left(-i\\frac{ω}{c}\\,Δ\\,d\\right) swap23.
```

`Δ = construct_Δ(k_par, construct_constitutive(ε, μ), a)` is built from the full 3×3 ε (with any
Euler rotation) and μ (the scalar fallback `μ_i·I`, or the layer's tensor `mu`).
Because Δ's eigenvalues are the mode wavevectors `q`, `exp(-i(ω/c)Δd)` equals the
eigenmode propagator `D·P(d)·D⁻¹` **without** diagonalizing, so this path needs no
eigenvalue sorting and is degeneracy-immune. `swap23` reorders between Berreman's
field vector `(Eₓ, Hᵧ, Eᵧ, −Hₓ)` and the dynamical-matrix basis; the `exp(-iωt)`
sign matches [`propagation_matrix`](@ref).

The matrix exponential is evaluated by scaling-and-squaring with degree-13 Padé
approximants (StaticArrays' `exp`).

Berreman, 1972, https://doi.org/10.1364/JOSA.62.000502
Mackay & Lakhtakia, 2020, https://doi.org/10.1007/978-3-031-02022-3
Higham, 2005, https://doi.org/10.1137/04061101X
"""
function layer_transfer_exp(layer, λ, k_par, ω, μ_i)
    ε = _layer_epsilon(layer, λ)
    μ = ismagnetic(layer) ? get_permeability(layer, λ) : permeability_tensor(μ_i, μ_i, μ_i)
    C = construct_constitutive(ε, μ)
    a = construct_a(k_par, C)
    Δ = construct_Δ(k_par, C, a)
    return _swap23 * exp(-im * (ω / c_0) * Δ * layer.thickness) * _swap23
end


"""
    construct_constitutive(ε, μ, ρ1, ρ2)

Construct the 6×6 constitutive matrix C = [ε ρ1; ρ2 μ] from the dielectric and permeability tensors.
"""
function construct_constitutive(ε, μ=Diagonal(ones(3)), ρ1=zeros(3, 3), ρ2=zeros(3, 3))
    return [ε ρ1; ρ2 μ]
end

# Specialized method for isotropic materials without magnetoelectric coupling (ρ₁ = ρ₂ = 0)
# This is the common case and avoids heap allocations from zeros() and hvcat
function construct_constitutive(ε::Diagonal{ComplexF64,SVector{3,ComplexF64}},
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
function construct_constitutive(ε::SMatrix{3,3,ComplexF64},
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

# Full 3×3 ε AND full 3×3 μ (e.g. gyrotropic/gyromagnetic media, no magnetoelectric coupling)
function construct_constitutive(ε::SMatrix{3,3,ComplexF64}, μ::SMatrix{3,3,ComplexF64})
    z = zero(ComplexF64)
    return @SMatrix [
        ε[1,1] ε[1,2] ε[1,3] z z z;
        ε[2,1] ε[2,2] ε[2,3] z z z;
        ε[3,1] ε[3,2] ε[3,3] z z z;
        z z z μ[1,1] μ[1,2] μ[1,3];
        z z z μ[2,1] μ[2,2] μ[2,3];
        z z z μ[3,1] μ[3,2] μ[3,3]
    ]
end

# Diagonal ε with a full μ tensor
construct_constitutive(ε::Diagonal{ComplexF64,SVector{3,ComplexF64}}, μ::SMatrix{3,3,ComplexF64}) =
    construct_constitutive(SMatrix{3,3,ComplexF64}(ε), μ)


"""
    construct_a(k_par, C)

Construct the elements of the intermediate 6×6 matrix ``a`` in terms of the
elements of matrix ``C`` (the 6×6 constitutive matrix holding the material dielectric and permeability tensors)
and propagation vector k_par. This is implemented as described in

Berreman, 1972, https://doi.org/10.1364/JOSA.62.000502
"""
function construct_a(k_par, C)
    d = C[3,3] * C[6,6] - C[3,6] * C[6,3]

    a31 = (C[6,1] * C[3,6] - C[3,1] * C[6,6]) / d
    a32 =((C[6,2] - k_par) * C[3,6] - C[3,2] * C[6,6]) / d
    a34 = (C[6,4] * C[3,6] -  C[3,4] * C[6,6]) / d
    a35 = (C[6,5] * C[3,6] - (C[3,5] + k_par) * C[6,6]) / d

    a61 = (C[6,3] * C[3,1] - C[3,3] * C[6,1]) / d
    a62 = (C[6,3] * C[3,2] - C[3,3] * (C[6,2] - k_par)) / d
    a64 = (C[6,3] * C[3,4] - C[3,3] * C[6,4]) / d
    a65 = (C[6,3] * (C[3,5] + k_par) - C[3,3] * C[6,5]) / d
    
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
    construct_Δ(k_par, C, a)

Construct the reordered matrix Δ in terms of the elements of
the two matrices, C and a, and the in-plane reduced wavevector k_par = ``k_x / k_0``.
The matrix Δ is involved in the relation

```math
    \\frac{\\delta}{\\delta z}\\Psi = \\frac{i \\omega}{c}\\Delta \\Psi
```

and Δ is the reordered S matrix in Berreman's formulation.

Berreman, 1972, https://doi.org/10.1364/JOSA.62.000502
"""
function construct_Δ(k_par, C, a)

    Δ11 =  C[5,1] + (C[5,3] + k_par) * a[3,1] + C[5,6] * a[6,1]
    Δ12 =  C[5,5] + (C[5,3] + k_par) * a[3,5] + C[5,6] * a[6,5]
    Δ13 =  C[5,2] + (C[5,3] + k_par) * a[3,2] + C[5,6] * a[6,2]
    Δ14 = -C[5,4] - (C[5,3] + k_par) * a[3,4] - C[5,6] * a[6,4]

    Δ21 =  C[1,1] + C[1,3] * a[3,1] + C[1,6] * a[6,1]
    Δ22 =  C[1,5] + C[1,3] * a[3,5] + C[1,6] * a[6,5]
    Δ23 =  C[1,2] + C[1,3] * a[3,2] + C[1,6] * a[6,2]
    Δ24 = -C[1,4] - C[1,3] * a[3,4] - C[1,6] * a[6,4]

    Δ31 = -C[4,1] - C[4,3] * a[3,1] - C[4,6] * a[6,1]
    Δ32 = -C[4,5] - C[4,3] * a[3,5] - C[4,6] * a[6,5]
    Δ33 = -C[4,2] - C[4,3] * a[3,2] - C[4,6] * a[6,2]
    Δ34 =  C[4,4] + C[4,3] * a[3,4] + C[4,6] * a[6,4]

    Δ41 =  C[2,1] + C[2,3] * a[3,1] + (C[2,6] - k_par) * a[6,1]
    Δ42 =  C[2,5] + C[2,3] * a[3,5] + (C[2,6] - k_par) * a[6,5]
    Δ43 =  C[2,2] + C[2,3] * a[3,2] + (C[2,6] - k_par) * a[6,2]
    Δ44 = -C[2,4] - C[2,3] * a[3,4] - (C[2,6] - k_par) * a[6,4]

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

    # Determine whether eigenvalues are effectively real. The eigen decomposition
    # of a real-valued Δ can produce tiny imaginary parts (≈ε_mach) due to floating
    # point arithmetic. Using `isreal()` (exact zero check) misclassifies these as
    # complex and sorts by imaginary part, which fails for nearly-real eigenvalues.
    effectively_real = all(abs(imag(q)) < 1e-10 * max(abs(q), 1.0) for q in q_unsorted)

    kt = 1
    kr = 1
    if effectively_real
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
    calculate_E_modes(k_par, q, ε, μ)

The 4 x 3 matrix E_modes contains vector components that belong
to the electric field calculated such that singularities are identified and removed.

`q[1]` and `q[2]` are forward-traveling modes and
`q[3]` and `q[4]` are backward-traveling modes.

Each row `v` of `E_modes` is normalized to unit length using the Hermitian norm
``‖v‖ = √(Σ|vᵢ|²)`` as prescribed by the 2019 erratum (Passler &
Paarmann, JOSAB 36, 3246). The Hermitian norm is the correct choice
because each row represents an electric field polarization direction: its
physical "length" is the field amplitude, which is ``√(E⋅E*) = √(Σ|Eᵢ|²)``.
The bilinear form ``√(Σ vᵢ²)`` has no physical meaning for complex fields
and can vanish for nonzero vectors (e.g. circular polarization),
making it unsuitable as a normalization factor.

This is based on the work in
Xu, et al., 2000, https://doi.org/10.1103/PhysRevB.61.1740
"""
function calculate_E_modes(k_par, q, ε, μ)

    E_modes = @MMatrix zeros(ComplexF64, 4, 3)
    E_modes[1,1] = 1
    E_modes[2,2] = 1
    E_modes[4,2] = 1
    E_modes[3,1] = -1

    # Common denominator term - can become singular at specific angles
    denom_33 = μ * ε[3,3] - k_par^2
    singular_33 = abs(denom_33) < eps(Float64)

    if isapprox(q[1], q[2])
        E_modes[1,2] = 0
        E_modes[2,1] = 0
    else
        E_modes[1,2] = (μ * ε[2,3] * (μ * ε[3,1] + k_par * q[1]) - μ * ε[2,1] * denom_33) / (denom_33 * (μ * ε[2,2] - k_par^2 - q[1]^2) - μ^2 * ε[2,3] * ε[3,2])
        E_modes[2,1] = (μ * ε[3,2] * (μ * ε[1,3] + k_par * q[2]) - μ * ε[1,2] * denom_33) / (denom_33 * (μ * ε[1,1] - q[2]^2) - (μ * ε[1,3] + k_par * q[2]) * (μ * ε[3,1] + k_par * q[2]))
    end

    # E_modes[i,3] components use denom_33 directly - set to zero if singular
    if singular_33
        E_modes[1,3] = 0
        E_modes[2,3] = 0
    else
        E_modes[1,3] = (-μ * ε[3,1] - k_par * q[1] - μ * ε[3,2] * E_modes[1,2]) / denom_33
        E_modes[2,3] = (-(μ * ε[3,1] + k_par * q[2]) * E_modes[2,1] - μ * ε[3,2]) / denom_33
    end

    if isapprox(q[3], q[4])
        E_modes[3,2] = 0.0
        E_modes[4,1] = 0.0
    else
        E_modes[3,2] = (μ * ε[2,1] * denom_33 - μ * ε[2,3] * (μ * ε[3,1] + k_par * q[3])) / (denom_33 * (μ * ε[2,2] - k_par^2 - q[3]^2) - μ^2 * ε[2,3] * ε[3,2])
        E_modes[4,1] = (μ * ε[3,2] * (μ * ε[1,3] + k_par * q[4]) - μ * ε[1,2] * denom_33) / (denom_33 * (μ * ε[1,1] - q[4]^2) - (μ * ε[1,3] + k_par * q[4]) * (μ * ε[3,1] + k_par * q[4]))
    end

    # E_modes[i,3] components for backward modes - set to zero if singular
    if singular_33
        E_modes[3,3] = 0
        E_modes[4,3] = 0
    else
        E_modes[3,3] = (μ * ε[3,1] + k_par * q[3] - μ * ε[3,2] * E_modes[3,2]) / denom_33
        E_modes[4,3] = (-(μ * ε[3,1] + k_par * q[4]) * E_modes[4,1] - μ * ε[3,2] ) / denom_33
    end

    # Normalize E_modes using the Hermitian norm: |v| = √(Σ|vᵢ|²).
    #
    # Each E_modes row is a polarization eigenvector of the electric field. The
    # physical amplitude of a complex field E is √(E⋅E*) = √(Σ|Eᵢ|²),
    # not √(Σ Eᵢ²). The bilinear form √(Σvᵢ²) conflates phase information
    # with amplitude — for example, v = [1, i, 0] (circular polarization)
    # gives Σvᵢ² = 0, a division-by-zero singularity, even though the field
    # has unit amplitude. The Hermitian norm is positive-definite on nonzero
    # vectors by construction, so it is always safe as a normalization factor.
    #
    # For real E_modes (lossless isotropic media) both norms coincide, so the
    # distinction only matters for complex E_modes (absorbing or anisotropic media).
    for i in 1:4
        v = SVector(E_modes[i,1], E_modes[i,2], E_modes[i,3])
        Z = norm(v)
        if abs(Z) > eps(Float64)
            E_modes[i,1] /= Z
            E_modes[i,2] /= Z
            E_modes[i,3] /= Z
        end
    end

    return SMatrix(E_modes)
end


# k̄ × v cross-product matrix for k̄ = (k_par, 0, q)
_kcross(k_par, q) = SMatrix{3,3,ComplexF64}(0, q, 0,  -q, 0, k_par,  0, -k_par, 0)


"""
    calculate_E_modes_tensor(k_par, q, ε, μ)

E-field eigenvectors for a layer with a full 3×3 permeability tensor `μ`, as the
null space of `W(q) = K μ⁻¹ K + ε` for each mode `q`. Degenerate pairs (qᵢ ≈ qⱼ,
e.g. isotropic media) share a 2-D null space; both orthogonal null vectors are
assigned, avoiding the singular dynamical matrix a naive per-mode null space
produces. Rows are unit-Hermitian-normalized, matching [`calculate_E_modes`].
"""
function calculate_E_modes_tensor(k_par, q, ε, μ; rtol=1e-7)
    μinv = inv(SMatrix{3,3,ComplexF64}(μ)); εm = SMatrix{3,3,ComplexF64}(ε)
    E_modes = @MMatrix zeros(ComplexF64, 4, 3); assigned = MVector{4,Bool}(false, false, false, false)
    for m in 1:4
        assigned[m] && continue
        K = _kcross(k_par, q[m]); F = svd(Matrix(K*μinv*K + εm))
        partner = 0
        for n in (m+1):4
            if !assigned[n] && abs(q[n]-q[m]) < rtol*max(abs(q[m]), 1.0)
                partner = n; break
            end
        end
        v1 = SVector{3,ComplexF64}(F.V[:,3]...); E_modes[m,:] = v1 ./ norm(v1); assigned[m] = true
        if partner != 0
            v2 = SVector{3,ComplexF64}(F.V[:,2]...); E_modes[partner,:] = v2 ./ norm(v2); assigned[partner] = true
        end
    end
    return SMatrix(E_modes)
end


"""
    dynamical_matrix(k_par, q, E_modes, μ::AbstractMatrix)

Dynamical matrix for a tensor μ. Rows are `(Eₓ, Eᵧ, Hᵧ, −Hₓ)` per mode with
`H = μ⁻¹ (k̄ × E)`. Reduces to the scalar-μ method when `μ = μ·I`.
"""
function dynamical_matrix(k_par, q, E_modes, μ::AbstractMatrix)
    μinv = inv(SMatrix{3,3,ComplexF64}(μ))
    cols = ntuple(4) do m
        E = SVector{3,ComplexF64}(E_modes[m,1], E_modes[m,2], E_modes[m,3])
        H = μinv * (_kcross(k_par, q[m]) * E)
        SVector{4,ComplexF64}(E[1], E[2], H[2], -H[1])
    end
    return hcat(cols...)::SMatrix{4,4,ComplexF64}
end


"""
    dynamical_matrix(k_par, q, E_modes, μ)

The dynamical matrix relating two layers at the interface
where matrix ``A_i`` for layer ``i`` relates the field ``E_i`` to
the field in the previous layer ``i - 1`` via

```math
    A_{i-1}E_{i-1} = A_{i}E_{i}
```

Xu, et al., 2000, https://doi.org/10.1103/PhysRevB.61.1740
"""
function dynamical_matrix(k_par, q, E_modes, μ)
    A_3 = (E_modes[:, 1] .* q .- k_par * E_modes[:, 3]) ./ μ
    A_4 = E_modes[:, 2] .* q ./ μ

    return @SMatrix [
        E_modes[1, 1] E_modes[2, 1] E_modes[3, 1] E_modes[4, 1];
        E_modes[1, 2] E_modes[2, 2] E_modes[3, 2] E_modes[4, 2];
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
