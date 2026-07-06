# Linear (p,s) -> circular (L,R) change-of-basis matrices for the exp(-iωt) time
# convention. Columns of C are ordered (L, R): e_L = (x̂ + iŷ)/√2, e_R = (x̂ - iŷ)/√2,
# so [Ex; Ey] = C [E_L; E_R]. The sign of i is fixed by exp(-iωt) (do NOT flip it
# without also changing the package time convention).
const _C_CIRC = SMatrix{2,2,ComplexF64}(1, im, 1, -im) / sqrt(2)      # [1 1; im -im]
const _C_CIRC_INV = SMatrix{2,2,ComplexF64}(1, 1, -im, im) / sqrt(2)  # [1 -im; 1 im]

# Transmission-null tolerance: below this the transmitted amplitude weight
# |t_c[L,j]|² + |t_c[R,j]|² is taken as zero, so circular transmittance is
# defined as 0 (avoids 0/0 at a perfect mirror).
const _CIRC_T_TOL = 1e-12

# Transform a 2×2 linear Jones matrix ([out,in] ordered p,s) into the circular
# (L,R) basis: J_circ = C⁻¹ J_lin C. Returns a 2×2 ([out,in] ordered L,R).
_jones_to_circular(J_lin) = _C_CIRC_INV * J_lin * _C_CIRC

# Assemble circular-basis R/T from the linear complex Jones coefficients
# r = (rpp,rps,rss,rsp), t = (tpp,tps,tsp,tss) and the Poynting flux factors:
# f_out1/f_out2 are the unit-amplitude z-fluxes of the two substrate eigenmodes
# (p-like / s-like), f_in_p/f_in_s those of the incident modes.
#
# Reflectance is |r_circ|^2. Transmittance for circular input j is computed by
# mapping the input Jones vector through t_lin to the substrate (p,s) amplitude
# pair `a`, whose per-mode flux `f_out1|a₁|² + f_out2|a₂|²` is the EXACT total
# transmitted power (distinct lossless eigenmodes carry no cross z-flux). That
# total is then split between the L/R outputs in proportion to |t_circ|². For an
# isotropic substrate (f_out1 = f_out2) the split itself is exact, because the
# circular unit vectors are orthonormal combinations of degenerate modes; for an
# anisotropic substrate the two eigenmodes are not circular, so only the per-input
# TOTAL is exact and the L/R split is the natural |t_circ|²-weighted estimate
# (see `transfer` docstring).
function _circular_result(r, t, f_out1, f_out2, f_in_p, f_in_s)
    # Package coefficients are indexed in-first (rps = p input → s output), so
    # the [out,in] Jones matrices place the p→s conversion in row 2, column 1.
    # (Using the transposed orientation is invisible for isotropic stacks —
    # their Jones matrices are diagonal — but breaks flux unitarity, and hence
    # the circular energy budget, for polarization-converting stacks.)
    r_lin = SMatrix{2,2,ComplexF64}(r[1], r[2], r[4], r[3])  # [rpp rsp; rps rss]
    t_lin = SMatrix{2,2,ComplexF64}(t[1], t[2], t[3], t[4])  # [tpp tsp; tps tss]
    r_c = _jones_to_circular(r_lin)
    t_c = _jones_to_circular(t_lin)
    Rc = abs2.(r_c)

    # Unit-amplitude circular incident flux: |χ_p|² f_in_p + |χ_s|² f_in_s with
    # χ = (1, ±i)/√2 (incident p/s modes carry no cross z-flux).
    f_in = (f_in_p + f_in_s) / 2

    Tc = MMatrix{2,2,Float64}(undef)
    for j in 1:2
        a = t_lin * _C_CIRC[:, j]                 # substrate (p,s) amplitudes for input j (1=L, 2=R)
        w = abs2(t_c[1, j]) + abs2(t_c[2, j])     # = |a₁|² + |a₂|² (C is unitary)
        if w > _CIRC_T_TOL
            T_tot = (f_out1 * abs2(a[1]) + f_out2 * abs2(a[2])) / f_in
            Tc[1, j] = T_tot * abs2(t_c[1, j]) / w
            Tc[2, j] = T_tot * abs2(t_c[2, j]) / w
        else
            Tc[1, j] = 0.0
            Tc[2, j] = 0.0
        end
    end
    # circular matrices are [out,in] ordered (L,R): [1,1]=LL [1,2]=LR [2,1]=RL [2,2]=RR
    return CircularTransferResult(Tc[2, 2], Tc[1, 1], Tc[2, 1], Tc[1, 2],
                                  Rc[2, 2], Rc[1, 1], Rc[2, 1], Rc[1, 2])
end


"""
    calculate_tr(M_sys)

Read the reflection and transmission coefficients off the full-stack transfer
matrix `M_sys`.

`M_sys` is the ``4\\times4`` transfer matrix for the whole stack, reordered into the
p/s block convention (the two p modes first, then the two s modes) so the
reflection and transmission Jones blocks can be read straight off its elements.
The function returns the co- and cross-polarized amplitude coefficients
(`rpp, rss, rps, rsp`, `tpp, tss, tps, tsp`) and their squared magnitudes.

Reflectance is ``R = |r|^2`` directly. The returned ``|t|^2`` values (all four,
co- and cross-polarized) are *not* physical transmittances — the transmitted
wave is in a different medium — so every transmittance reported by `transfer`
(`Tpp`, `Tss`, `Tps`, `Tsp`) comes from the per-mode Poynting-vector
calculation (`calculate_tr(S::Poynting)`) instead.

The ``2\\times2`` transfer-matrix relations originate with Yeh (1979); the
generalized ``4\\times4`` coefficient formulas follow Passler & Paarmann (2017).

References:
- Yeh, 1979, https://doi.org/10.1364/JOSA.69.000742
- Passler & Paarmann, 2017, https://doi.org/10.1364/JOSAB.34.002128
"""
function calculate_tr(M_sys)

    d = M_sys[1,1] * M_sys[3,3] - M_sys[1,3] * M_sys[3,1]

    rpp = (M_sys[2,1] * M_sys[3,3] - M_sys[2,3] * M_sys[3,1]) / d
    rss = (M_sys[1,1] * M_sys[4,3] - M_sys[4,1] * M_sys[1,3]) / d
    rps = (M_sys[4,1] * M_sys[3,3] - M_sys[4,3] * M_sys[3,1]) / d
    rsp = (M_sys[1,1] * M_sys[2,3] - M_sys[2,1] * M_sys[1,3]) / d

    tpp =  M_sys[3,3] / d
    tss =  M_sys[1,1] / d
    tps = -M_sys[3,1] / d
    tsp = -M_sys[1,3] / d

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
