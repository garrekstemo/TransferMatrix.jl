# Independent reference oracles for validating anisotropic-μ physics (issue #91).
#
# These are deliberately INDEPENDENT of the package's eigenvector/γ path:
#
#   * `transfer_exp`  — matrix-exponential propagator. Each interior layer's
#     transfer matrix is `Λ · exp(-i(ω/c)Δ·d) · Λ`, built from the Berreman Δ
#     matrix (which already supports a full 3×3 μ tensor). Needs no γ formulas
#     and is degeneracy-immune. This is the reference the #91 implementation
#     must reproduce, and the recommended implementation route.
#
#   * `abeles_RT`     — closed-form 2×2 Abelès characteristic-matrix solution
#     for ISOTROPIC magnetic multilayers (external analytic ground truth).
#
#   * `transfer_candidate` (gamma_nullspace_robust) — an eigenvector-γ
#     implementation using the generalized dynamical matrix `H = μ⁻¹(k̄×E)`,
#     with degeneracy-robust null-space handling. Cross-checks the exp route.
#
# Conventions match the package exactly (exp(-iωt); field vector ordering
# Ψ = (Ex, Hy, Ey, -Hx); the Λ₁₃₂₄ reorder).
module MuReference

using TransferMatrix, LinearAlgebra, StaticArrays
const TM = TransferMatrix
const C0 = TM.c_0
const Λ1324 = @SMatrix [1.0 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

export GLayer, iso, glay, vac, transfer_exp, transfer_candidate, exp_rt,
       glayer_from_package, gamma_nullspace_robust, abeles_RT, media_to_glayers,
       μ_iso, μ_diag, μ_polder, μ_rot

# A general layer: full ε(λ) and μ(λ) tensors + thickness.
struct GLayer
    εf::Function
    μf::Function
    d::Float64
end
iso(n, d=0.0) = GLayer(λ -> SMatrix{3,3,ComplexF64}(n^2,0,0, 0,n^2,0, 0,0,n^2),
                       λ -> SMatrix{3,3,ComplexF64}(I), Float64(d))
glay(εt, μt, d) = GLayer(λ -> εt, λ -> μt, Float64(d))
vac() = iso(1.0)

μ_iso(m)        = SMatrix{3,3,ComplexF64}(m,0,0, 0,m,0, 0,0,m)
μ_diag(a,b,c)   = SMatrix{3,3,ComplexF64}(a,0,0, 0,b,0, 0,0,c)
μ_polder(m,κ)   = SMatrix{3,3,ComplexF64}(m,-im*κ,0, im*κ,m,0, 0,0,1)   # bias ∥ z
μ_rot(a,b,c,φ,θ,ψ) = (R = TM.euler_rotation_matrix(φ,θ,ψ);
                      SMatrix{3,3,ComplexF64}(R*Matrix(μ_diag(a,b,c))*R'))

# k̄ × v  cross-product matrix for k̄ = (ξ, 0, q)
Kmat(ξ, q) = SMatrix{3,3,ComplexF64}(0, q, 0,  -q, 0, ξ,  0, -ξ, 0)

# Boundary dynamical matrix for an ISOTROPIC half-space (μ scalar·I), using the
# package's validated scalar-μ γ path (degeneracy-safe for isotropic media).
function boundaryD(L::GLayer, λ, ξ)
    ε = SMatrix{3,3,ComplexF64}(L.εf(λ)); μ = SMatrix{3,3,ComplexF64}(L.μf(λ))
    μs = μ[1,1]
    M = TM.construct_constitutive(Diagonal(ε), Diagonal(μ))
    a = TM.construct_a(ξ, M); Δ = TM.construct_Δ(ξ, M, a)
    q, _ = TM.calculate_q(Δ, a); q = ComplexF64.(q)
    γ = TM.calculate_E_modes(ξ, q, Diagonal(ε), μs)
    D = TM.dynamical_matrix(ξ, q, γ, μs)
    return Matrix(D), SMatrix{4,3,ComplexF64}(γ), q
end

# Interior-layer transfer (D-convention) via the matrix exponential of Δ.
function layerT_exp(L::GLayer, λ, ξ, ω)
    ε = SMatrix{3,3,ComplexF64}(L.εf(λ)); μ = SMatrix{3,3,ComplexF64}(L.μf(λ))
    M = TM.construct_constitutive(ε, μ); a = TM.construct_a(ξ, M); Δ = TM.construct_Δ(ξ, M, a)
    return Matrix(Λ1324) * exp(-im * (ω / C0) * Matrix(Δ) * L.d) * Matrix(Λ1324)
end

# Tensor-μ dynamical matrix: rows (Ex, Ey, Hy, -Hx) with H = μ⁻¹(k̄×E).
function dynmat_tensor(μ, ξ, q, γ)
    μinv = inv(SMatrix{3,3,ComplexF64}(μ))
    cols = ntuple(4) do m
        E = SVector{3,ComplexF64}(γ[m,1], γ[m,2], γ[m,3])
        H = μinv * (Kmat(ξ, q[m]) * E)
        SVector{4,ComplexF64}(E[1], E[2], H[2], -H[1])
    end
    return hcat(cols...)::SMatrix{4,4,ComplexF64}
end

# Degeneracy-robust E-eigenvectors as the null space of W(q) = K μ⁻¹ K + ε.
function gamma_nullspace_robust(ε, μ, ξ, q; rtol=1e-7)
    μinv = inv(SMatrix{3,3,ComplexF64}(μ)); εm = SMatrix{3,3,ComplexF64}(ε)
    γ = zeros(ComplexF64, 4, 3); assigned = falses(4)
    for m in 1:4
        assigned[m] && continue
        K = Kmat(ξ, q[m]); F = svd(Matrix(K*μinv*K + εm))
        partner = 0
        for n in (m+1):4
            if !assigned[n] && abs(q[n]-q[m]) < rtol*max(abs(q[m]),1.0)
                partner = n; break
            end
        end
        if partner == 0
            γ[m,:] = F.V[:,3]/norm(F.V[:,3]); assigned[m] = true
        else
            γ[m,:] = F.V[:,3]/norm(F.V[:,3]); γ[partner,:] = F.V[:,2]/norm(F.V[:,2])
            assigned[m] = true; assigned[partner] = true
        end
    end
    return SMatrix{4,3,ComplexF64}(γ)
end

# Full stack solve via the matrix-exponential propagator. Returns R/T (+Γ, r/t).
function transfer_exp(stack::Vector{GLayer}, λ; θ=0.0)
    N = length(stack)
    ε_amb = SMatrix{3,3,ComplexF64}(stack[1].εf(λ))
    ξ = sqrt(ε_amb[1,1]) * sin(θ); ω = 2π * C0 / λ
    D1, γ1, q1 = boundaryD(stack[1], λ, ξ)
    DN, γN, qN = boundaryD(stack[N], λ, ξ)
    Γ = inv(D1)
    for i in 2:N-1
        Γ = Γ * layerT_exp(stack[i], λ, ξ, ω)
    end
    Γ = Matrix(Λ1324) * (Γ * DN) * Matrix(Λ1324)
    Γs = SMatrix{4,4,ComplexF64}(Γ)
    r, R, t, T = TM.calculate_tr(Γs)
    S = TM.poynting(ξ, q1, qN, γ1, γN, t, r)
    Tpp, Tss, _, _, Tps, Tsp = TM.calculate_tr(S)
    return (; Γ=Γs, r=r, t=t,
            Rpp=R[1], Rss=R[2], Rsp=R[3], Rps=R[4],
            Tpp=Tpp, Tss=Tss, Tps=Tps, Tsp=Tsp)
end

# Full stack solve using a candidate γ for interior layers (boundaries: package).
function transfer_candidate(stack::Vector{GLayer}, λ; θ=0.0, gammafn=gamma_nullspace_robust)
    N = length(stack)
    ε_amb = SMatrix{3,3,ComplexF64}(stack[1].εf(λ))
    ξ = sqrt(ε_amb[1,1]) * sin(θ); ω = 2π * C0 / λ
    function layerD(L)
        ε = SMatrix{3,3,ComplexF64}(L.εf(λ)); μ = SMatrix{3,3,ComplexF64}(L.μf(λ))
        M = TM.construct_constitutive(ε, μ); a = TM.construct_a(ξ, M); Δ = TM.construct_Δ(ξ, M, a)
        q, _ = TM.calculate_q(Δ, a); q = ComplexF64.(q)
        D = dynmat_tensor(μ, ξ, q, gammafn(ε, μ, ξ, q))
        return Matrix(D), TM.propagation_matrix(ω, q)
    end
    D1, _, _ = boundaryD(stack[1], λ, ξ)
    Γ = inv(D1)
    for i in 2:N-1
        Di, Pi = layerD(stack[i])
        Γ = Γ * (Di * Matrix(Pi(stack[i].d)) * inv(Di))
    end
    DN, _, _ = boundaryD(stack[N], λ, ξ)
    Γ = Matrix(Λ1324) * (Γ * DN) * Matrix(Λ1324)
    r, _, t, _ = TM.calculate_tr(SMatrix{4,4,ComplexF64}(Γ))
    return r, t
end

exp_rt(stack, λ; θ=0.0) = (g = transfer_exp(stack, λ; θ=θ); (g.r, g.t))

# Build a GLayer stack mirroring a package Layer stack with a scalar μ everywhere.
function glayer_from_package(layer::Layer, μscalar::Real)
    function εf(λ)
        nx, ny, nz = TM.get_refractive_indices(layer, λ)
        εd = TM.dielectric_tensor(TM.dielectric_constant(nx), TM.dielectric_constant(ny), TM.dielectric_constant(nz))
        φ, θe, ψ = TM.get_euler_angles(layer)
        if φ != 0.0 || θe != 0.0 || ψ != 0.0
            return SMatrix{3,3,ComplexF64}(TM.rotate_dielectric_tensor(εd, TM.euler_rotation_matrix(φ, θe, ψ)))
        end
        return SMatrix{3,3,ComplexF64}(εd)
    end
    return GLayer(εf, λ -> SMatrix{3,3,ComplexF64}(μscalar * I), Float64(layer.thickness))
end

# ---- Analytic Abelès oracle for ISOTROPIC magnetic media (vacuum ambient) ----
abeles_layer(η, δ) = [cos(δ) -im*sin(δ)/η; -im*η*sin(δ) cos(δ)]   # exp(-iωt)
function abeles_RT(media, λ, θ, pol)
    ξ = sin(θ); k0 = 2π/λ
    qof(ε,μ) = sqrt(ε*μ - ξ^2); ηof(ε,μ) = pol === :s ? qof(ε,μ)/μ : ε/qof(ε,μ)
    ε0,μ0,_ = media[1]; εs,μs,_ = media[end]; η0 = ηof(ε0,μ0); ηsub = ηof(εs,μs)
    M = ComplexF64[1 0; 0 1]
    for i in 2:length(media)-1
        ε,μ,d = media[i]; M = M * abeles_layer(ηof(ε,μ), k0*qof(ε,μ)*d)
    end
    m11,m12,m21,m22 = M[1,1],M[1,2],M[2,1],M[2,2]
    num = η0*m11 + η0*ηsub*m12 - m21 - ηsub*m22
    den = η0*m11 + η0*ηsub*m12 + m21 + ηsub*m22
    r = num/den; t = 2η0/den
    return abs2(r), real(ηsub)/real(η0)*abs2(t)
end
media_to_glayers(media) = [GLayer(let v=ε; λ->SMatrix{3,3,ComplexF64}(v,0,0,0,v,0,0,0,v) end,
                                   let v=μ; λ->SMatrix{3,3,ComplexF64}(v,0,0,0,v,0,0,0,v) end,
                                   Float64(d)) for (ε,μ,d) in media]

end # module
