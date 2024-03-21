using Test, TransferMatrix
# Test types.jl

@testset "Layer" begin
    l = TransferMatrix.Layer(10e-6, 1.000273, 0.0)
    @test l.t ≈ 0.00001
    @test l.n == 1.000273
    @test l.κ == 0.0
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
#     @test s.layers[1].n == fill(1.0, length(λs))
#     @test length(s.layers[2].n) == length(λs)
#     @test length(s.layers[2].κ) == length(λs)
# end