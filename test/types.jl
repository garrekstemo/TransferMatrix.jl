# Test types.jl

@testset "Layer" begin
    au = RefractiveMaterial("main", "Au", "Rakic-LD")
    l_au = TransferMatrix.Layer(au, 10e-6)
    
    @test l_au.thickness ≈ 0.00001
    @test dispersion(l_au.material, 2.0) == 0.777149018404908
    @test extinction(l_au.material, 2.0) == 12.597128834355829
    
    sio2 = RefractiveMaterial("main", "SiO2", "Malitson")
    @test_throws DomainError TransferMatrix.Layer(sio2, -1.0)
end