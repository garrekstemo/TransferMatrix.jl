# Incoherent (thick) layers by intensity-matrix chaining.
#
# A layer much thicker than the coherence length of the light (a mm-thick window or
# substrate measured with an FTIR, whose etalon fringes are finer than the
# instrument resolution) adds intensities rather than amplitudes. The stack is cut
# at each incoherent layer into coherent blocks; every block is solved with the 2×2
# path for its forward and backward power coefficients (Rf, Tf, Rb, Tb), and the
# blocks are chained through the incoherent layers' single-pass power attenuation
# A = exp(−2 k₀ Im(q) d), with q = n cos θ the layer's forward wavevector.
#
# Katsidis & Siapkas, Appl. Opt. 41, 3978 (2002), https://doi.org/10.1364/AO.41.003978
# write the chain as a product of 2×2 intensity matrices
#   (1/Tf) [1  −Rb; Rf  Tf·Tb − Rf·Rb]  and  diag(1/A, A),
# with T = 1/P₁₁ and R = P₂₁/P₁₁ for the product P. Here the same product is
# accumulated as the equivalent pairwise recursion
#   Tf = Tf₁ A Tf₂ / (1 − Rb₁ A² Rf₂),   Rf = Rf₁ + Tf₁ Tb₁ A² Rf₂ / (1 − Rb₁ A² Rf₂),
# (and likewise backward), which never forms 1/A and so stays finite for strongly
# absorbing thick layers.

_incoherent_indices(i::Integer) = [Int(i)]
_incoherent_indices(v) = sort!(unique!(Int[Int(i) for i in v]))

function _validate_incoherent(idx, layers, method, basis)
    N = length(layers)
    for i in idx
        (2 ≤ i ≤ N - 1) || throw(ArgumentError(
            "incoherent layer index $i out of range; must be an interior layer, 2 ≤ i ≤ $(N - 1)"))
        _isotropic_boundary(layers[i]) || throw(ArgumentError(
            "incoherent layer $i must be isotropic and non-magnetic"))
    end
    _decoupled_structure(layers) || throw(ArgumentError(
        "incoherent layers need an isotropic, non-magnetic ambient and substrate"))
    method === :auto || method === :decoupled || throw(ArgumentError(
        "incoherent layers are computed with the 2×2 path; use method=:auto or :decoupled, got :$method"))
    basis === :linear || throw(ArgumentError(
        "incoherent layers support basis=:linear only, got :$basis"))
    return nothing
end

# Forward and backward power coefficients of the coherent block layers[a:b], whose
# end layers are incoherent media (or the ambient/substrate), for p and s.
function _block_powers(::Type{C}, λ, layers, a, b, k_par, μ, sd, ε_a, ε_b) where {C}
    res = _block_matrices_2x2(C, λ, layers, a, b, k_par, μ, sd, ε_a, ε_b)
    res isa _Widen && return res
    res === nothing && throw(ArgumentError(
        "the coherent block between layers $a and $b couples s and p; incoherent " *
        "layers need every coherent block to be decoupled (see `isdecoupled`)"))
    Mp, Ms, med_a, med_b = res
    # forward-mode z-fluxes on either side (backward modes carry the same magnitude)
    return _powers(Mp, med_a.Sp_f[3], med_b.Sp_f[3]), _powers(Ms, med_a.Ss_f[3], med_b.Ss_f[3])
end

# (Rf, Tf, Rb, Tb) from a 2×2 mode-amplitude transfer matrix [a] = M [b] and the
# unit-amplitude fluxes of the forward modes on either side.
@inline function _powers(M, f_a, f_b)
    tf = 1 / M[1,1]
    rf = M[2,1] * tf
    rb = -M[1,2] * tf
    tb = (M[1,1] * M[2,2] - M[1,2] * M[2,1]) * tf
    return (abs2(rf), abs2(tf) * f_b / f_a, abs2(rb), abs2(tb) * f_a / f_b)
end

# Chain (Rf₁, Tf₁, Rb₁, Tb₁), single-pass attenuation A, then (Rf₂, Tf₂, Rb₂, Tb₂).
@inline function _chain_powers(P1, A, P2)
    Rf1, Tf1, Rb1, Tb1 = P1
    Rf2, Tf2, Rb2, Tb2 = P2
    A2 = A * A
    den = 1 - Rb1 * A2 * Rf2
    return (Rf1 + Tf1 * Tb1 * A2 * Rf2 / den,
            Tf1 * A * Tf2 / den,
            Rb2 + Tb2 * Tf2 * A2 * Rb1 / den,
            Tb2 * A * Tb1 / den)
end

function _transfer_incoherent(::Type{C}, λ, layers, θ, μ, sd, idx, validate) where {C}
    N = length(layers)
    k0 = (2π * c_0 / λ) / c_0
    bounds = [1; idx; N]
    ok_a, ε_a = _boundary_epsilon(C, layers[1], λ)
    ok_b, ε_b = _boundary_epsilon(C, layers[bounds[2]], λ)
    (ok_a && ok_b) || return _WIDEN
    k_par = √ε_a * sin(θ)
    res = _block_powers(C, λ, layers, 1, bounds[2], k_par, μ, sd, ε_a, ε_b)
    res isa _Widen && return res
    Pp, Ps = res
    for j in 2:length(bounds) - 1
        m = bounds[j]
        ε_m = ε_b
        ok_b, ε_b = _boundary_epsilon(C, layers[bounds[j + 1]], λ)
        ok_b || return _WIDEN
        q_m = _forward_branch(sqrt(complex(μ * ε_m - k_par^2)))
        A = exp(-2 * k0 * imag(q_m) * layers[m].thickness)
        res = _block_powers(C, λ, layers, m, bounds[j + 1], k_par, μ, sd, ε_m, ε_b)
        res isa _Widen && return res
        Qp, Qs = res
        Pp = _chain_powers(Pp, A, Qp)
        Ps = _chain_powers(Ps, A, Qs)
    end
    Rpp, Tpp = Pp[1], Pp[2]
    Rss, Tss = Ps[1], Ps[2]
    z = zero(Tpp)
    if validate
        _validate_physics(λ, layers, Tpp, Tss, Rpp, Rss; sheets=sd)
    end
    return TransferResult(promote(Tpp, Tss, z, z, Rpp, Rss, z, z)...)
end
