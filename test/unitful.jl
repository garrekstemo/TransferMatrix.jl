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
        @test TransferMatrix._to_radians(0.5) === 0.5
        @test TransferMatrix._to_um.([0.1, 0.2]) == [0.1, 0.2]
        @test TransferMatrix._to_wavelength_um.([1.0, 2.0]) == [1.0, 2.0]
        @test TransferMatrix._to_radians.([0.0, 0.3]) == [0.0, 0.3]
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

    @testset "angle of incidence accepts angle units" begin
        @test TransferMatrix._to_radians(45u"°") ≈ deg2rad(45)
        @test TransferMatrix._to_radians((π / 4)u"rad") ≈ π / 4

        ref = transfer(1.55, layers; θ=deg2rad(45))
        @test tr_equal(transfer(1.55, layers; θ=45u"°"), ref)
        @test tr_equal(transfer(1.55, layers; θ=(π / 4)u"rad"),
                       transfer(1.55, layers; θ=π / 4))
    end

    @testset "sweep_angle accepts a vector of angle units" begin
        λs = [1.4, 1.55, 1.7]
        θs_deg = [0u"°", 30u"°", 60u"°"]
        θs_rad = deg2rad.([0.0, 30.0, 60.0])
        @test sweep_angle(λs, θs_deg, layers).Rpp ≈ sweep_angle(λs, θs_rad, layers).Rpp
    end

    @testset "field profiles accept angle units" begin
        Edeg = efield(1.55, layers; θ=30u"°")
        Erad = efield(1.55, layers; θ=deg2rad(30))
        @test Edeg.p ≈ Erad.p
        @test Edeg.s ≈ Erad.s
    end

    @testset "sweeps accept unit-bearing vectors" begin
        λs = [1.4, 1.55, 1.7]
        θs = [0.0, 0.3]
        @test sweep_angle(λs .* u"μm", θs, layers).Rpp ≈ sweep_angle(λs, θs, layers).Rpp

        ts = [0.1, 0.2, 0.3]
        ru = sweep_thickness(λs .* u"μm", ts .* u"μm", layers, 2)
        rp = sweep_thickness(λs, ts, layers, 2)
        @test ru.Rpp ≈ rp.Rpp
        @test ru.Tss ≈ rp.Tss
    end

    @testset "tabulated-data Layer accepts unit-bearing thickness" begin
        λs_data = [1.0, 1.5, 2.0]
        ns_data = [1.5, 1.5, 1.5]
        ks_data = [0.0, 0.0, 0.0]
        layer = Layer(λs_data, ns_data, ks_data, 100u"nm")
        @test layer.thickness ≈ 0.1
        @test layer.thickness isa Float64
        @test !isanisotropic(layer)
        stack = [Layer(n_air, 0.1), layer, Layer(n_sub, 0.5)]
        @test transfer(1.5u"μm", stack) isa TransferResult
    end

    @testset "field profiles accept unit-bearing λ and dz" begin
        Eu = efield(1.55u"μm", layers; dz=1u"nm")
        Ep = efield(1.55, layers; dz=0.001)
        @test Eu.z ≈ Ep.z
        @test Eu.p ≈ Ep.p
        @test Eu.s ≈ Ep.s

        Hu = hfield(1.55u"μm", layers; dz=1u"nm")
        Hp = hfield(1.55, layers; dz=0.001)
        @test Hu.z ≈ Hp.z
        @test Hu.p ≈ Hp.p
    end

end
