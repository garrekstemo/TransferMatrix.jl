# Linear (p,s) -> circular (L,R) change-of-basis matrices for the exp(-iωt) time
# convention. Columns of C are ordered (L, R): e_L = (x̂ + iŷ)/√2, e_R = (x̂ - iŷ)/√2,
# so [Ex; Ey] = C [E_L; E_R]. The sign of i is fixed by exp(-iωt) (do NOT flip it
# without also changing the package time convention).
const _C_CIRC = SMatrix{2,2,ComplexF64}(1, im, 1, -im) / sqrt(2)      # [1 1; im -im]
const _C_CIRC_INV = SMatrix{2,2,ComplexF64}(1, 1, -im, im) / sqrt(2)  # [1 -im; 1 im]

# Transmission-null tolerance: below this the transmitted Poynting flux is taken
# as zero, so circular transmittance is defined as 0 (avoids 0/0 at a perfect mirror).
const _CIRC_T_TOL = 1e-12

# Transform a 2×2 linear Jones matrix ([out,in] ordered p,s) into the circular
# (L,R) basis: J_circ = C⁻¹ J_lin C. Returns a 2×2 ([out,in] ordered L,R).
_jones_to_circular(J_lin) = _C_CIRC_INV * J_lin * _C_CIRC

# Assemble circular-basis R/T from the linear complex Jones coefficients
# r = (rpp,rps,rss,rsp), t = (tpp,tps,tsp,tss) and the Poynting diagonal
# transmittances Tpp, Tss. Reflectance is |r_circ|^2; transmittance is N·|t_circ|^2
# with the single polarization-independent Poynting scalar N (exact for an isotropic
# substrate; see `transfer` docstring for caveats).
function _circular_result(r, t, Tpp, Tss)
    r_lin = SMatrix{2,2,ComplexF64}(r[1], r[4], r[2], r[3])  # [rpp rps; rsp rss]
    t_lin = SMatrix{2,2,ComplexF64}(t[1], t[3], t[2], t[4])  # [tpp tps; tsp tss]
    r_c = _jones_to_circular(r_lin)
    t_c = _jones_to_circular(t_lin)
    Rc = abs2.(r_c)
    denom = abs2(t[1]) + abs2(t[4])  # |tpp|^2 + |tss|^2
    N = denom > _CIRC_T_TOL ? (Tpp + Tss) / denom : zero(denom)
    Tc = N .* abs2.(t_c)
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

Reflectance is ``R = |r|^2`` directly. The returned ``|t|^2`` is *not* the
physical transmittance — the transmitted wave is in a different medium — so the
diagonal `Tpp`/`Tss` reported by `transfer` come from the Poynting-vector
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
