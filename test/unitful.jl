using Test
using TransferMatrix
using Unitful

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

    tr_equal(a, b) = all(f -> getfield(a, f) ≈ getfield(b, f),
                         (:Tpp, :Tss, :Tps, :Tsp, :Rpp, :Rss, :Rps, :Rsp))

    @testset "Layer accepts unit-bearing thickness" begin
        layer = Layer(n_film, 100u"nm")
        @test layer.thickness ≈ 0.1
        @test layer.thickness isa Float64
        aniso = Layer(n_film, n_film, n_sub, 250u"nm")
        @test aniso.thickness ≈ 0.25
        @test aniso.thickness isa Float64
    end

    @testset "transfer accepts unit-bearing wavelength (length)" begin
        ref = transfer(1.55, layers)
        @test tr_equal(transfer(1.55u"μm", layers), ref)
        @test tr_equal(transfer(1550u"nm", layers), ref)
    end

    @testset "transfer accepts spectral wavelength inputs" begin
        ref = transfer(1.55, layers)
        λ = 1.55u"μm"
        @test tr_equal(transfer(1 / λ, layers), ref)                       # wavenumber ν̃ = 1/λ
        @test tr_equal(transfer(Unitful.c0 / λ, layers), ref)              # frequency  λ = c/f
        @test tr_equal(transfer(Unitful.h * Unitful.c0 / λ, layers), ref)  # energy     λ = hc/E
    end

    @testset "non-spectral units are rejected" begin
        @test_throws ArgumentError transfer(5u"kg", layers)
    end

end
