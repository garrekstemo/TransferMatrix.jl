# Test types.jl

@testset "Layer" begin
    n_au = RefractiveMaterial("main", "Au", "Rakic-LD")
    au = TransferMatrix.Layer(n_au, 10e-6)
    
    @test au.thickness ≈ 0.00001
    @test real(au.dispersion(2.0)) == 0.777149018404908
    @test imag(au.dispersion(2.0)) == 12.597128834355829
    
    n_sio2 = RefractiveMaterial("main", "SiO2", "Malitson")
    
    @test isa(au, TransferMatrix.Layer)
    @test_throws DomainError TransferMatrix.Layer(n_sio2, -1.0)

    λs = [1.0, 2.0, 3.0]
    dispersion = [1.5, 2.0, 2.5]
    extinction = [0.5, 1.0, 1.5]
    thickness = 1.0
    layer = TransferMatrix.Layer(λs, dispersion, extinction, thickness)
    
    @test isa(layer, TransferMatrix.Layer)
    @test layer.thickness == thickness
    @test layer.dispersion(1.0) == 1.5 + 0.5im
    @test layer.dispersion(2.0) == 2.0 + 1.0im
    @test layer.dispersion(3.0) == 2.5 + 1.5im

    zero_layer = TransferMatrix.Layer(λs, dispersion, extinction, 0.0)
    @test zero_layer.thickness == 0.0
end
