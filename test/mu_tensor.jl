# Anisotropic-μ physics validation (issue #91) and regression guards for the
# anisotropic mode-sorting fix in `evaluate_birefringence`.
#
# Until a μ-tensor API lands, the *reference* implementations live in
# `mu_reference.jl` (an independent matrix-exponential propagator + an analytic
# Abelès oracle + a degeneracy-robust eigenvector-γ route). These tests:
#   1. regression-test the SHIPPING package against the independent oracles, and
#   2. lock in the validated μ-tensor physics the #91 implementation must match.
# When #91 is implemented, point the `transfer_exp`/`transfer_candidate` calls at
# the public API to turn these into direct API tests.

include("mu_reference.jl")
using .MuReference

const MU_TOL = 1e-9      # machine-precision agreement (results observed ~1e-15)
const E_TOL  = 1e-7      # energy / duality tolerance

@testset "μ-tensor physics (issue #91)" begin

    # ---------------------------------------------------------------------
    # (1) REGRESSION: axis-aligned biaxial/uniaxial must be finite AND energy-
    # conserving at oblique incidence. Guards the evaluate_birefringence fix
    # (previously NaN at θ>0 due to a degenerate Poynting p/s ratio).
    # ---------------------------------------------------------------------
    @testset "axis-aligned anisotropy: finite + energy-conserving (sorting fix)" begin
        amb = Layer(λ -> 1.0, 1.0); sub = Layer(λ -> 1.0, 1.0)
        biax = Layer(λ -> 2.0, λ -> 2.5, λ -> 3.0, 0.5)        # axis-aligned biaxial
        uni  = Layer(λ -> 1.6, λ -> 1.6, λ -> 1.9, 0.4)        # axis-aligned uniaxial
        for L in (biax, uni), θ in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
            r = transfer(1.0, [amb, L, sub]; θ=θ)
            for v in (r.Rpp, r.Rss, r.Rps, r.Rsp, r.Tpp, r.Tss, r.Tps, r.Tsp)
                @test isfinite(v)
            end
            @test isapprox(r.Rpp + r.Rps + r.Tpp, 1.0; atol=E_TOL)  # p-input energy
            @test isapprox(r.Rss + r.Rsp + r.Tss, 1.0; atol=E_TOL)  # s-input energy
        end
    end

    # ---------------------------------------------------------------------
    # (2) The shipping package's scalar-μ path matches the independent
    # matrix-exponential reference to machine precision (incl. magnetic μ≠1,
    # rotated anisotropy, absorption, oblique incidence).
    # ---------------------------------------------------------------------
    @testset "package scalar-μ == matrix-exponential reference" begin
        amb = Layer(λ -> 1.0, 1.0); sub = Layer(λ -> 1.5, 1.0)
        glass = Layer(λ -> 1.5, 0.3)
        uni  = Layer(λ -> 1.6, λ -> 1.6, λ -> 1.8, 0.4; euler=(0.4, 0.6, 0.2))
        biax = Layer(λ -> 1.5, λ -> 1.7, λ -> 1.9, 0.35; euler=(0.5, 0.7, 0.3))
        absL = Layer(λ -> 1.5 + 0.05im, 0.3)
        for θ in (0.0, 0.3, 0.6, 0.9), μ in (1.0, 1.4, 2.0)
            for ls in ([amb, glass, sub], [amb, uni, sub], [amb, biax, sub],
                       [amb, absL, sub], [amb, glass, uni, biax, sub])
                ref = transfer(1.0, ls; θ=θ, μ=μ)
                g = transfer_exp([glayer_from_package(L, μ) for L in ls], 1.0; θ=θ)
                @test isapprox(ref.Rpp, g.Rpp; atol=MU_TOL)
                @test isapprox(ref.Rss, g.Rss; atol=MU_TOL)
                @test isapprox(ref.Tpp, g.Tpp; atol=MU_TOL)
                @test isapprox(ref.Tss, g.Tss; atol=MU_TOL)
            end
        end
    end

    # ---------------------------------------------------------------------
    # (3) Matrix-exponential reference == analytic magnetic Abelès oracle
    # (isotropic magnetic media; external closed-form ground truth).
    # ---------------------------------------------------------------------
    @testset "matrix-exp reference == analytic magnetic Abelès" begin
        slab  = [(1.0+0im,1.0+0im,0.0),(4.0+0im,2.5+0im,0.35),(1.0+0im,1.0+0im,0.0)]
        lossy = [(1.0+0im,1.0+0im,0.0),(4.0+0.2im,2.5+0.1im,0.3),(1.5+0im,1.0+0im,0.0)]
        multi = [(1.0+0im,1.0+0im,0.0),(4.0+0im,2.0+0im,0.2),(2.25+0im,3.0+0im,0.25),
                 (6.0+0im,1.5+0im,0.15),(2.0+0im,1.0+0im,0.0)]
        for θ in (0.0, 0.25, 0.6, 1.0), pol in (:s, :p), med in (slab, lossy, multi)
            Ra, Ta = abeles_RT(med, 1.0, θ, pol)
            g = transfer_exp(media_to_glayers(med), 1.0; θ=θ)
            Rg = pol === :s ? g.Rss : g.Rpp
            Tg = pol === :s ? g.Tss : g.Tpp
            @test isapprox(Ra, Rg; atol=MU_TOL)
            @test isapprox(Ta, Tg; atol=MU_TOL)   # μ_sub == μ_amb for these stacks
        end
    end

    # ---------------------------------------------------------------------
    # (4) μ-tensor: eigenvector-γ (H = μ⁻¹(k̄×E)) == matrix-exp reference, for
    # diagonal, rotated, gyromagnetic (Polder) and lossy-gyromagnetic μ.
    # Confirms the dynamical-matrix generalization the #91 implementation needs.
    # ---------------------------------------------------------------------
    @testset "eigenvector-γ (H=μ⁻¹ k×E) == matrix-exp, for tensor μ" begin
        εg = μ_iso(2.25)
        cases = [μ_diag(2,2,3), μ_diag(1.5,2,3), μ_rot(1.5,2,3,0.4,0.6,0.3),
                 μ_polder(2.0,0.6),
                 SMatrix{3,3,ComplexF64}(2+0.1im,-im*0.5,0, im*0.5,2+0.1im,0, 0,0,1)]
        for θ in (0.0, 0.2, 0.5, 0.9), μt in cases
            st = [vac(), glay(εg, μt, 0.3), vac()]
            re, te = exp_rt(st, 1.0; θ=θ)
            rc, tc = transfer_candidate(st, 1.0; θ=θ)
            @test maximum(abs, rc .- re) < MU_TOL
            @test maximum(abs, tc .- te) < MU_TOL
        end
    end

    # ---------------------------------------------------------------------
    # (5) Electric-magnetic DUALITY (#86 ↔ #91): a purely magnetic slab
    # (ε=I, μ=tensor) via the exp reference equals the dual purely-dielectric
    # slab (ε=tensor, μ=1) via the PACKAGE's full-ε path, under s↔p swap.
    # This also exercises the package's internal full-3×3-ε machinery.
    # ---------------------------------------------------------------------
    @testset "electric-magnetic duality (#86 ↔ #91)" begin
        for θ in (0.0, 0.2, 0.5, 0.9), (a,b,c,φ,ϑ,ψ) in
                ((1.5,2.0,3.0,0.4,0.6,0.3), (1.2,1.8,2.5,0.7,0.5,0.9))
            orig = [vac(), glay(μ_iso(1.0), μ_rot(a,b,c,φ,ϑ,ψ), 0.3), vac()]
            go = transfer_exp(orig, 1.0; θ=θ)
            dual = [Layer(λ->1.0,1.0),
                    Layer(λ->sqrt(a), λ->sqrt(b), λ->sqrt(c), 0.3; euler=(φ,ϑ,ψ)),
                    Layer(λ->1.0,1.0)]
            gd = transfer(1.0, dual; θ=θ, μ=1.0)
            @test isapprox(go.Rpp, gd.Rss; atol=E_TOL)
            @test isapprox(go.Rss, gd.Rpp; atol=E_TOL)
            @test isapprox(go.Rps, gd.Rsp; atol=E_TOL)
            @test isapprox(go.Tpp, gd.Tss; atol=E_TOL)
        end
    end

    # ---------------------------------------------------------------------
    # (6) Energy conservation for lossless μ-tensor stacks (vacuum/vacuum):
    # Rpp+Rps+Tpp+Tps = 1 and Rss+Rsp+Tss+Tsp = 1 (each transmittance is the
    # per-output-mode Poynting flux of its own channel). Includes gyromagnetic
    # (Hermitian Polder) μ, which converts p↔s in transmission.
    # ---------------------------------------------------------------------
    @testset "energy conservation, lossless μ-tensor" begin
        for θ in (0.0, 0.2, 0.5, 0.9),
            μt in (μ_diag(1.5,2.5,3.5), μ_rot(1.5,2.5,3.5,0.4,0.6,0.3),
                   μ_polder(2.0,0.6), μ_polder(2.5,1.2))
            g = transfer_exp([vac(), glay(μ_iso(2.25), μt, 0.4), vac()], 1.0; θ=θ)
            @test isapprox(g.Rpp + g.Rps + g.Tpp + g.Tps, 1.0; atol=E_TOL)
            @test isapprox(g.Rss + g.Rsp + g.Tss + g.Tsp, 1.0; atol=E_TOL)
        end
    end

    # ---------------------------------------------------------------------
    # (7) Gyromagnetic non-reciprocity (Faraday): the transmission Jones matrix
    # is antisymmetric (t_ps = -t_sp, vs +t_sp for reciprocal media), the
    # Onsager relation t_sp(κ) = t_ps(-κ) holds, and cross-pol vanishes at κ=0.
    # ---------------------------------------------------------------------
    @testset "gyromagnetic non-reciprocity (Faraday)" begin
        slab(κ) = [vac(), glay(μ_iso(2.25), μ_polder(2.0, κ), 0.4), vac()]
        for θ in (0.0, 0.4)
            _, t = exp_rt(slab(0.6), 1.0; θ=θ)
            @test isapprox(t[2], -t[3]; atol=MU_TOL)            # t_ps = -t_sp
            @test abs(t[2]) > 1e-3                               # genuinely nonzero
        end
        _, tp = exp_rt(slab(0.6), 1.0; θ=0.4)
        _, tm = exp_rt(slab(-0.6), 1.0; θ=0.4)
        @test isapprox(tp[3], tm[2]; atol=MU_TOL)               # Onsager t_sp(κ)=t_ps(-κ)
        _, t0 = exp_rt(slab(0.0), 1.0; θ=0.4)
        @test abs(t0[2]) < MU_TOL && abs(t0[3]) < MU_TOL        # κ=0 ⇒ no cross-pol
    end
end

@testset "public transfer honors per-layer μ == exp reference" begin
    host(μt) = [Layer(λ->1.0, 1.0), Layer(λ->1.5, 0.3; mu=μt), Layer(λ->1.0, 1.0)]
    refstack(μt) = [iso(1.0), glay(μ_iso(2.25), μt, 0.3), iso(1.0)]
    cases = [μ_diag(2,2,3), μ_rot(1.5,2,3,0.4,0.6,0.3), μ_polder(2.0,0.6)]
    for θ in (0.0, 0.3, 0.7), μt in cases
        r = transfer(1.0, host(μt); θ=θ)
        g = transfer_exp(refstack(μt), 1.0; θ=θ)
        @test isapprox(r.Rpp, g.Rpp; atol=1e-9)
        @test isapprox(r.Rss, g.Rss; atol=1e-9)
        @test isapprox(r.Tpp, g.Tpp; atol=1e-9)
        @test isapprox(r.Rps, g.Rps; atol=1e-9)
        @test isapprox(r.Rsp, g.Rsp; atol=1e-9)
    end
end

@testset "sweep_thickness preserves per-layer μ" begin
    μt = μ_polder(2.0, 0.6)
    layers = [Layer(λ->1.0, 1.0), Layer(λ->1.5, 0.3; mu=μt), Layer(λ->1.0, 1.0)]
    ts = [0.2, 0.3, 0.4]
    λs = [0.9, 1.0, 1.1]
    sw = sweep_thickness(λs, ts, layers, 2)
    for (ti, t) in enumerate(ts)
        layers_t = [Layer(λ->1.0, 1.0), Layer(λ->1.5, t; mu=μt), Layer(λ->1.0, 1.0)]
        for (ji, λ) in enumerate(λs)
            r = transfer(λ, layers_t)
            @test isapprox(sw.Rpp[ti, ji], r.Rpp; atol=1e-12)
            @test isapprox(sw.Tpp[ti, ji], r.Tpp; atol=1e-12)
            @test isapprox(sw.Rps[ti, ji], r.Rps; atol=1e-12)
        end
    end
end

@testset "hfield finite + continuous for a magnetic layer" begin
    layers = [Layer(λ->1.0, 1.0), Layer(λ->1.5, 0.3; mu=gyrotropic_tensor(2.0, 0.6)), Layer(λ->1.0, 1.0)]
    H = hfield(1.0, layers; θ=0.3, dz=0.01)
    @test all(isfinite, H.p) && all(isfinite, H.s)
end
