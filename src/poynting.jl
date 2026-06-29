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
    poynting(k_par, q_in, q_out, γ_in, γ_out, t_coefs, r_coefs[, μ_in, μ_out])

Calculate the Poynting vector from wavevectors ``q``,
components of the electric field γ, and transmission
and reflection coefficients.

Transmitted Poynting vectors use substrate wavevectors (`q_out`), while
reflected Poynting vectors use incident-medium wavevectors (`k_in[3,:]`,
`k_in[4,:]`), since reflected waves propagate in the incident medium.

The optional `μ_in` and `μ_out` arguments are the 3×3 permeability tensors of
the ambient (incident) and substrate media respectively (default: identity, i.e.
non-magnetic). They are used to form ``H = μ^{-1}(k̄ \\times E)`` when computing
the Poynting flux, so that energy conservation holds for a magnetic substrate.
A magnetic ambient is only correctly handled at normal incidence; at oblique
incidence the conserved in-plane wavevector k_par is still computed from the
ambient permittivity alone (see issue #71), so results may be inaccurate.

!!! note "Transmittance vs reflectance"
    This function computes Poynting vectors for both transmitted and reflected
    waves, but only the **transmitted** Poynting vectors are used for the final
    output. Reflectance is computed as ``R = |r|^2`` from the transfer matrix
    coefficients — see [`transfer`](@ref) for the rationale.
"""
function poynting(k_par, q_in, q_out, γ_in, γ_out, t_coefs, r_coefs,
                  μ_in::AbstractMatrix  = SMatrix{3,3,ComplexF64}(I),
                  μ_out::AbstractMatrix = SMatrix{3,3,ComplexF64}(I))

    μin_inv  = inv(SMatrix{3,3,ComplexF64}(μ_in))
    μout_inv = inv(SMatrix{3,3,ComplexF64}(μ_out))

    # create the wavevector in the first layer
    k_in = @MMatrix zeros(ComplexF64, 4, 3)
    k_in[:, 1] .= k_par

    for (i, q_i) in enumerate(q_in)
        k_in[i, 3] = q_i
    end

    k_in ./= c_0
    k_in = SMatrix(k_in)

    E_forward_in_p =  γ_in[1, :]  # p-polarized incident electric field
    E_forward_in_s =  γ_in[2, :]  # s-polarized incident electric field
    # E_backward_in_p = γ_in[3, :]
    # E_backward_in_s = γ_in[4, :]

    # Each transmitted/reflected field is a superposition of the two
    # substrate (resp. incident) eigenmodes, which carry the substrate mode-1
    # field γ_out[1] and mode-2 field γ_out[2].
    E_out_p1 = t_coefs[1] * γ_out[1, :]
    E_out_p2 = t_coefs[2] * γ_out[2, :]
    E_out_s1 = t_coefs[3] * γ_out[1, :]
    E_out_s2 = t_coefs[4] * γ_out[2, :]

    E_ref_p = r_coefs[1] * γ_in[3, :] + r_coefs[2] * γ_in[4, :]
    E_ref_s = r_coefs[3] * γ_in[3, :] + r_coefs[4] * γ_in[4, :]

    S_in_p = real(0.5 * E_forward_in_p × conj(μin_inv * (k_in[1, :] × E_forward_in_p)))
    S_in_s = real(0.5 * E_forward_in_s × conj(μin_inv * (k_in[2, :] × E_forward_in_s)))

    k_out = @MMatrix zeros(ComplexF64, 4, 3)
    k_out[:, 1] .= k_par

    for (i, q_i) in enumerate(q_out)
        k_out[i, 3] = q_i
    end

    k_out ./= c_0
    k_out = SMatrix(k_out)

    # Transmitted power is the sum of the two substrate eigenmodes' Poynting
    # vectors, EACH evaluated with its OWN wavevector (k_out modes 1, 2). The
    # two modes carry different wavevectors when the substrate is anisotropic,
    # and the cross-interference between distinct propagating modes transports
    # no net z-flux, so the per-mode sum is the physical transmitted flux.
    # (For an isotropic substrate k_out[1]=k_out[2] and this reduces to the
    # single-wavevector form, since the p/s cross term carries no z-power.)
    S_out_p = real(0.5 * E_out_p1 × conj(μout_inv * (k_out[1, :] × E_out_p1))) +
              real(0.5 * E_out_p2 × conj(μout_inv * (k_out[2, :] × E_out_p2)))
    S_out_s = real(0.5 * E_out_s1 × conj(μout_inv * (k_out[1, :] × E_out_s1))) +
              real(0.5 * E_out_s2 × conj(μout_inv * (k_out[2, :] × E_out_s2)))

    # Reflected waves: use incident-medium wavevectors (k_in modes 3,4),
    # not substrate wavevectors, because reflected light propagates backward
    # through the incident medium with its own q-values.
    S_refl_p = real(0.5 * E_ref_p × conj(μin_inv * (k_in[3, :] × E_ref_p)))
    S_refl_s = real(0.5 * E_ref_s × conj(μin_inv * (k_in[4, :] × E_ref_s)))

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
# Sort one mode pair (a transmitted or reflected pair) so the p-like mode is
# first and the s-like mode second, matching the fixed polarization references
# in `calculate_γ`.
#
# The Poynting ratio C = |Sx|²/(|Sx|²+|Sy|²) distinguishes p from s only when
# the modes actually differ along x vs y. For an axis-aligned (diagonal-ε)
# crystal BOTH eigenmodes have Sy = 0, so C ≈ 1 for both and the ratio is
# degenerate — it cannot tell p from s. In that case fall back to the
# electric-field ratio |Ex|²/(|Ex|²+|Ey|²), which separates the p-mode
# (Ex ≠ 0) from the s-mode (Ey ≠ 0). `isapprox(NaN, NaN) = false`, so a 0/0
# Poynting ratio (both Sx and Sy zero) is also routed to the E-field fallback
# via the `!isfinite` checks.
function sort_polarization_pair!(modes, Ψ, S)
    Cp = abs_ratio(S[1, modes[1]], S[2, modes[1]])
    Cs = abs_ratio(S[1, modes[2]], S[2, modes[2]])
    if isapprox(Cp, Cs) || !isfinite(Cp) || !isfinite(Cs)
        Cp = abs_ratio(Ψ[1, modes[1]], Ψ[3, modes[1]])
        Cs = abs_ratio(Ψ[1, modes[2]], Ψ[3, modes[2]])
    end
    if Cs > Cp
        reverse!(modes)
    end
    return modes
end

function evaluate_birefringence(Ψ, S, t_modes, r_modes)

    # Sort each pair (transmitted, reflected) so the p-like mode is first
    # (slot 1 / slot 3) and the s-like mode second (slot 2 / slot 4). This
    # ordering is assumed by the fixed polarization references in
    # `calculate_γ` (γ[1,1]=1, γ[2,2]=1, γ[3,1]=-1, γ[4,2]=1); getting it
    # wrong assigns each mode the wrong polarization, which both makes the
    # cross-polarization denominators in `calculate_γ` vanish (0/0 → NaN) and
    # corrupts the r/t coefficients (energy is not conserved).
    sort_polarization_pair!(t_modes, Ψ, S)
    sort_polarization_pair!(r_modes, Ψ, S)

    return t_modes, r_modes
end


"""
Ratio of the absolution squares of two components
used to evaluate if a material is birefringent.
"""
abs_ratio(a, b) = abs2(a) / (abs2(a) + abs2(b))


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
