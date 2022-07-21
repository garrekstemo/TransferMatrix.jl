using Test
using TransferMatrix

# Test dataio.jl

@testset "read_refractive" begin
    aufile = abspath("../refractive_index_data/Au_n_0.191-20.9um_Ciesielski.csv")
    caf2file = abspath("../refractive_index_data/CaF2_n_0.23-9.7um_Malitson.csv")
    au = TransferMatrix.read_refractive(aufile, "Au", 1.0, 1e6)
    caf2 = TransferMatrix.read_refractive(caf2file, "CaF2", 0.0, 1e6)

    @test au.thickness == 1.0
    @test au.material == "Au"
    @test length(au.λ) == 333
    @test au.λ[1] ≈ 0.19077e-6

    @test caf2.thickness == 0.0
    @test caf2.κ == fill(0.0, length(caf2.λ))
end