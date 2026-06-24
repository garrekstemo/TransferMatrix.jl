# Matrix-exponential propagator (issue #92): the public :exp backend, validated
# against the eigenmode path (:eig), the independent oracles in mu_reference.jl,
# energy conservation, and the interior-degeneracy robustness win.

if !isdefined(@__MODULE__, :MuReference)
    include("mu_reference.jl")
end
using .MuReference

@testset "matrix-exponential propagator (issue #92)" begin

    @testset "layer_transfer_exp == eig interior transfer (D-basis)" begin
        amb = Layer(λ -> 1.0, 1.0)
        cases = [
            ("isotropic",    Layer(λ -> 1.5, 0.3),                                      1.0),
            ("rot-uniaxial", Layer(λ -> 1.6, λ -> 1.6, λ -> 1.8, 0.4; euler=(0.4, 0.6, 0.2)), 1.0),
            ("scalar mu",    Layer(λ -> 1.5, 0.3),                                      2.0),
            ("tensor mu",    Layer(λ -> 1.5, 0.3; mu=gyrotropic_tensor(2.0, 0.6)),      1.0),
        ]
        for (_, L, μi) in cases, θ in (0.0, 0.6)
            λ = 1.0
            ω = 2π * TransferMatrix.c_0 / λ
            nx, _, _ = TransferMatrix.get_refractive_indices(amb, λ)
            ξ = sqrt(TransferMatrix.dielectric_constant(nx)) * sin(θ)
            D, P, _, _ = TransferMatrix.layer_matrices(L, λ, ξ, μi)
            Dm = SMatrix{4,4,ComplexF64}(D)
            T_eig = Dm * SMatrix{4,4,ComplexF64}(P(L.thickness)) * inv(Dm)
            T_exp = TransferMatrix.layer_transfer_exp(L, λ, ξ, ω, μi)
            @test maximum(abs, T_exp - T_eig) < 1e-12
        end
    end

    @testset "transfer :exp == :eig (no sheets)" begin
        amb = Layer(λ -> 1.0, 1.0); sub = Layer(λ -> 1.5, 1.0)
        glass = Layer(λ -> 1.5, 0.3)
        uni  = Layer(λ -> 1.6, λ -> 1.6, λ -> 1.8, 0.4; euler=(0.4, 0.6, 0.2))
        biax = Layer(λ -> 1.5, λ -> 1.7, λ -> 1.9, 0.35; euler=(0.5, 0.7, 0.3))
        absL = Layer(λ -> 1.5 + 0.05im, 0.3)
        magL = Layer(λ -> 1.5, 0.3; mu=gyrotropic_tensor(2.0, 0.6))
        stacks = [[amb, glass, sub], [amb, uni, sub], [amb, biax, sub],
                  [amb, absL, sub], [amb, glass, uni, biax, sub], [amb, magL, sub]]
        for ls in stacks, θ in (0.0, 0.3, 0.6, 0.9), μ in (1.0, 2.0)
            re = transfer(1.0, ls; θ=θ, μ=μ, method=:eig)
            rx = transfer(1.0, ls; θ=θ, μ=μ, method=:exp)
            for f in (:Tpp, :Tss, :Tps, :Tsp, :Rpp, :Rss, :Rps, :Rsp)
                @test isapprox(getfield(re, f), getfield(rx, f); atol=1e-12)
            end
        end
    end

    @testset "invalid method throws ArgumentError" begin
        ls = [Layer(λ -> 1.0, 1.0), Layer(λ -> 1.5, 1.0)]
        @test_throws ArgumentError transfer(1.0, ls; method=:bogus)
    end

    @testset "transfer :exp == :eig with conductive sheets" begin
        amb = Layer(λ -> 1.0, 1.0); sub = Layer(λ -> 1.5, 1.0)
        film = Layer(λ -> 1.6, 0.3)
        graphene = Sheet(6.08e-5 + 1.0e-5im)            # isotropic σ (S)
        tmdc = Sheet(; xx=1.0e-4 + 2.0e-5im, yy=8.0e-5 + 1.0e-5im)
        ls = [amb, film, sub]
        for θ in (0.0, 0.4, 0.8)
            for sh in (Dict(1 => graphene), Dict(2 => tmdc), Dict(1 => graphene, 2 => tmdc))
                re = transfer(1.0, ls; θ=θ, sheets=sh, method=:eig)
                rx = transfer(1.0, ls; θ=θ, sheets=sh, method=:exp)
                for f in (:Tpp, :Tss, :Tps, :Tsp, :Rpp, :Rss, :Rps, :Rsp)
                    @test isapprox(getfield(re, f), getfield(rx, f); atol=1e-12)
                end
            end
        end
    end

    @testset "sweeps honor method (:exp == :eig == per-call transfer)" begin
        amb = Layer(λ -> 1.0, 1.0)
        film = Layer(λ -> 1.6, λ -> 1.6, λ -> 1.8, 0.4; euler=(0.2, 0.5, 0.1))
        sub = Layer(λ -> 1.5, 1.0)
        layers = [amb, film, sub]
        λs = [0.9, 1.0, 1.1]; θs = [0.0, 0.4, 0.8]
        se = sweep_angle(λs, θs, layers; method=:eig, threads=false)
        sx = sweep_angle(λs, θs, layers; method=:exp, threads=false)
        for f in (:Tpp, :Rpp, :Rps, :Tss), i in eachindex(θs), j in eachindex(λs)
            @test isapprox(getfield(se, f)[i, j], getfield(sx, f)[i, j]; atol=1e-12)
            ref = transfer(λs[j], layers; θ=θs[i], method=:exp)
            @test isapprox(getfield(sx, f)[i, j], getfield(ref, f); atol=1e-12)
        end
        ts = [0.2, 0.4]
        twe = sweep_thickness(λs, ts, layers, 2; method=:eig, threads=false)
        twx = sweep_thickness(λs, ts, layers, 2; method=:exp, threads=false)
        for ti in eachindex(ts), j in eachindex(λs)
            @test isapprox(twe.Tpp[ti, j], twx.Tpp[ti, j]; atol=1e-12)
        end
    end

    @testset "circular basis under :exp == :eig" begin
        amb = Layer(λ -> 1.0, 1.0); sub = Layer(λ -> 1.4, 1.0)
        film = Layer(λ -> 1.5, 0.3; mu=gyrotropic_tensor(2.0, 0.6))
        for θ in (0.0, 0.5)
            ce = transfer(1.0, [amb, film, sub]; θ=θ, basis=:circular, method=:eig)
            cx = transfer(1.0, [amb, film, sub]; θ=θ, basis=:circular, method=:exp)
            for f in (:Trr, :Tll, :Trl, :Tlr, :Rrr, :Rll, :Rrl, :Rlr)
                @test isapprox(getfield(ce, f), getfield(cx, f); atol=1e-12)
            end
        end
    end

    @testset "energy conservation under :exp (lossless)" begin
        amb = Layer(λ -> 1.0, 1.0); sub = Layer(λ -> 1.0, 1.0)
        films = (Layer(λ -> 1.6, λ -> 1.6, λ -> 1.9, 0.4; euler=(π/6, π/4, 0)),  # out-of-plane tilt (#70)
                 Layer(λ -> 1.5, 0.3; mu=gyrotropic_tensor(2.0, 0.6)))           # gyromagnetic μ
        for film in films, θ in (0.0, 0.3, 0.6, 0.9)
            r = transfer(1.0, [amb, film, sub]; θ=θ, method=:exp)
            @test isapprox(r.Rpp + r.Rps + r.Tpp, 1.0; atol=1e-9)
            @test isapprox(r.Rss + r.Rsp + r.Tss, 1.0; atol=1e-9)
        end
    end

    @testset "interior mixed propagating/evanescent: :eig throws, :exp conserves" begin
        amb = Layer(λ -> 2.0, 1.0)                          # higher-index ambient
        interior = Layer(λ -> 1.4, λ -> 1.6, λ -> 1.7, 0.5) # axis-aligned, all n < 2
        sub = Layer(λ -> 2.5, 1.0)                          # clean (propagating) substrate
        stack = [amb, interior, sub]; θ = 0.97
        @test 1.6 < 2.0 * sin(θ) < 1.7                      # one transmitted mode evanescent
        @test_throws ArgumentError transfer(1.0, stack; θ=θ, method=:eig)
        r = transfer(1.0, stack; θ=θ, method=:exp)
        @test all(isfinite, (r.Rpp, r.Rps, r.Tpp, r.Rss, r.Rsp, r.Tss))
        @test isapprox(r.Rpp + r.Rps + r.Tpp, 1.0; atol=1e-9)
        @test isapprox(r.Rss + r.Rsp + r.Tss, 1.0; atol=1e-9)
    end

    @testset ":exp public API == MuReference.transfer_exp" begin
        amb = Layer(λ -> 1.0, 1.0); sub = Layer(λ -> 1.5, 1.0)
        uni = Layer(λ -> 1.6, λ -> 1.6, λ -> 1.8, 0.4; euler=(0.4, 0.6, 0.2))
        ls = [amb, uni, sub]
        for θ in (0.0, 0.6)
            rx = transfer(1.0, ls; θ=θ, method=:exp)
            g = transfer_exp([glayer_from_package(L, 1.0) for L in ls], 1.0; θ=θ)
            @test isapprox(rx.Rpp, g.Rpp; atol=1e-12)
            @test isapprox(rx.Tpp, g.Tpp; atol=1e-12)
            @test isapprox(rx.Rps, g.Rps; atol=1e-12)
        end
    end

    @testset "default method is :exp" begin
        # Stack where :eig throws but :exp works; the no-method default must NOT throw.
        amb = Layer(λ -> 2.0, 1.0)
        interior = Layer(λ -> 1.4, λ -> 1.6, λ -> 1.7, 0.5)
        sub = Layer(λ -> 2.5, 1.0)
        r = transfer(1.0, [amb, interior, sub]; θ=0.97)   # no method ⇒ default
        @test isapprox(r.Rpp + r.Rps + r.Tpp, 1.0; atol=1e-9)
        @test isapprox(r.Rss + r.Rsp + r.Tss, 1.0; atol=1e-9)
    end

end
