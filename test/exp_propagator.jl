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

end
