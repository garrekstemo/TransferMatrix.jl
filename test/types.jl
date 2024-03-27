using Test
using RefractiveIndex
using TransferMatrix

# Test types.jl

@testset "Layer" begin
    au = RefractiveMaterial("main", "Au", "Rakic-LD")
    l_au = TransferMatrix.Layer(au, 10e-6)
    
    @test l_au.thickness ≈ 0.00001
    @test dispersion(l_au.material, 2.0) == 0.777149018404908
    @test extinction(l_au.material, 2.0) == 12.597128834355829
    
    sio2 = RefractiveMaterial("main", "SiO2", "Malitson")
    @test_throws DomainError TransferMatrix.Layer(sio2, -1e-9)
end

# @testset "Structure" begin
    
#     l1 = TransferMatrix.Layer("Air", 0.0, [1e-9], [1.0], [0.0])

#     λ2 = [1.53, 1.54, 1.55, 1.56, 1.57] .* 1e-9
#     n2 = [1.443, 1.442, 1.441, 1.440, 1.399]
#     κ2 = [0.00002, 0.00002, 0.00002, 0.00002, 0.00002]
#     l2 = TransferMatrix.Layer("Glass", 0.0, λ2, n2, κ2)
    
#     θs = [0.0, 5.0, 10.0] .* π/180
#     λs = collect(1.531:0.001:1.566) .* 1e-9
#     s = TransferMatrix.Structure([l1, l2], λs, θs)

#     @test s.layers[1].λ == λs
#     @test s.layers[2].λ == λs
#     @test s.λ == λs
#     @test s.θ == θs
#     @test s.layers[1].n_re == fill(1.0, length(λs))
#     @test length(s.layers[2].n_re) == length(λs)
#     @test length(s.layers[2].n_im) == length(λs)
# end