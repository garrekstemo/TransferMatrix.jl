@testset "fresnel" begin
    n1 = 1.0  # air
    n2 = 1.5  # glass

    θ = 0.0
    Rs, Rp = fresnel(θ, n1, n2)

    @test Rs ≈ Rp  # at normal incidence, s and p are equal

    θ = π / 4
    Rs, Rp = fresnel(θ, n1, n2)

    @test Rp ≈ Rs^2  # relationship at 45 degrees

    θB = atan(n2 / n1)  # Brewster angle
    Rs, Rp = fresnel(θB, n1, n2)
    @test Rp ≈ 0.0 atol=1e-10  # p-polarized reflectance vanishes at Brewster's angle
    @test Rs > 0.0  # s-polarized reflectance remains nonzero

    # Test grazing incidence (θ → π/2)
    θ_grazing = π/2
    Rs, Rp = fresnel(θ_grazing, n1, n2)
    @test Rs == 1.0  # total reflection at grazing incidence
    @test Rp == 1.0

    # Test near-grazing incidence (approaches 1.0 but computed normally)
    θ_near_grazing = π/2 - 1e-10
    Rs, Rp = fresnel(θ_near_grazing, n1, n2)
    @test Rs ≈ 1.0 atol=1e-8
    @test Rp ≈ 1.0 atol=1e-8

    # Test total internal reflection (n1 > n2, θ > critical angle)
    n1_tir = 1.5  # glass
    n2_tir = 1.0  # air
    θ_critical = asin(n2_tir / n1_tir)  # critical angle ≈ 41.8°

    # Below critical angle: partial reflection
    θ_below = θ_critical - 0.1
    Rs, Rp = fresnel(θ_below, n1_tir, n2_tir)
    @test Rs < 1.0
    @test Rp < 1.0

    # Above critical angle: total internal reflection
    θ_above = θ_critical + 0.1
    Rs, Rp = fresnel(θ_above, n1_tir, n2_tir)
    @test Rs == 1.0
    @test Rp == 1.0

    # Well above critical angle
    θ_well_above = π/3  # 60°, well above ~41.8° critical angle
    Rs, Rp = fresnel(θ_well_above, n1_tir, n2_tir)
    @test Rs == 1.0
    @test Rp == 1.0
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
