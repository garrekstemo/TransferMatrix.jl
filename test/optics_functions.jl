@testset "fresnel" begin
    n1 = 1.0  # air
    n2 = 1.5  # glass

    θ = 0.0
    Rs, Rp = fresnel(θ, n1, n2)
    
    @test Rs ≈ Rp  # expected reflectance for s-polarized light

    θ = π / 4
    Rs, Rp = fresnel(θ, n1, n2)

    @test Rp ≈ Rs^2  # expected reflectance for s-polarized light
    # @test Rp ≈ 0.11  # expected reflectance for p-polarized light

    θ = π / 2  # 90 degrees
    Rs, Rp = fresnel(θ, n1, n2)

    @test Rs ≈ 1.0  # expected reflectance for s-polarized light at Brewster's angle
    @test Rp ≈ 1.0  # expected reflectance for p-polarized light at Brewster's angle
end

@testset "stopband" begin
    Δf = TransferMatrix.stopband(2.5, 1.5) * 630  # for 630 nm light
    @test isapprox(Δf, 202.685, atol = 1e-3)
end

@testset "dbr_reflectivity" begin
    no = 1.0  # air
    ns = 1.5  # glass
    n1 = 2.5  # TiO2
    n2 = 1.5  # SiO2

    r = dbr_reflectivity(no, ns, n1, n2, 3)
    @test isapprox(r, 0.883, atol = 1e-3)

    r = dbr_reflectivity(no, ns, n1, n2, 6)
    @test isapprox(r, 0.994, atol = 1e-3)
    
    r = dbr_reflectivity(no, ns, n1, n2, 9)
    @test isapprox(r, 0.999, atol = 1e-3)
end