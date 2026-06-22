using Test
using TransferMatrix

@testset "Unitful support" begin

    n_air = λ -> 1.0 + 0.0im
    n_film = λ -> 1.5 + 0.0im
    n_sub = λ -> 1.45 + 0.0im
    layers = [Layer(n_air, 0.1), Layer(n_film, 0.25), Layer(n_sub, 0.5)]

    @testset "core converter seams are no-ops for Real" begin
        @test TransferMatrix._to_um(0.1) === 0.1
        @test TransferMatrix._to_wavelength_um(1.55) === 1.55
        @test TransferMatrix._to_um.([0.1, 0.2]) == [0.1, 0.2]
        @test TransferMatrix._to_wavelength_um.([1.0, 2.0]) == [1.0, 2.0]
    end

    @testset "wiring leaves plain-number results unchanged" begin
        r = transfer(1.55, layers)
        @test r isa TransferResult
        @test 0.0 ≤ r.Rpp ≤ 1.0
        @test r.Rpp + r.Tpp ≈ 1.0 atol=1e-6
    end

end
