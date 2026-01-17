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
    @test Rs ≈ 1.0 atol=1e-10
    @test Rp ≈ 1.0 atol=1e-10

    # Well above critical angle
    θ_well_above = π/3  # 60°, well above ~41.8° critical angle
    Rs, Rp = fresnel(θ_well_above, n1_tir, n2_tir)
    @test Rs ≈ 1.0 atol=1e-10
    @test Rp ≈ 1.0 atol=1e-10
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

@testset "fresnel_coefficients" begin
    n1 = 1.0  # air
    n2 = 1.5  # glass

    # Normal incidence: check consistency with fresnel()
    θ = 0.0
    rs, rp, ts, tp = fresnel_coefficients(θ, n1, n2)
    Rs, Rp = fresnel(θ, n1, n2)

    @test abs2(rs) ≈ Rs
    @test abs2(rp) ≈ Rp
    @test abs(rs) ≈ abs(rp)  # s and p magnitudes equal at normal incidence (signs differ by convention)

    # Check transmission coefficients at normal incidence
    # t = 2n1 / (n1 + n2)
    expected_t = 2 * n1 / (n1 + n2)
    @test real(ts) ≈ expected_t
    @test real(tp) ≈ expected_t

    # Energy conservation at normal incidence: R + T = 1
    # T = (n2/n1) * |t|²  at normal incidence
    Ts = (n2 / n1) * abs2(ts)
    Tp = (n2 / n1) * abs2(tp)
    @test Rs + Ts ≈ 1.0 atol=1e-10
    @test Rp + Tp ≈ 1.0 atol=1e-10

    # Brewster angle: rp = 0
    θB = atan(n2 / n1)
    rs, rp, ts, tp = fresnel_coefficients(θB, n1, n2)
    @test abs(rp) ≈ 0.0 atol=1e-10
    @test abs(rs) > 0.0

    # Total internal reflection: |r| = 1 with phase shift
    n1_tir = 1.5
    n2_tir = 1.0
    θ_above_critical = π / 3  # 60°, above critical angle ~41.8°
    rs, rp, ts, tp = fresnel_coefficients(θ_above_critical, n1_tir, n2_tir)
    @test abs(rs) ≈ 1.0 atol=1e-10
    @test abs(rp) ≈ 1.0 atol=1e-10
    @test ts == 0.0  # no transmitted power
    @test tp == 0.0
end

@testset "airy" begin
    # Test 1: Single interface (d = 0) should match Fresnel
    n0 = 1.0
    nf = 1.5  # film
    ns = 1.5  # substrate same as film = effectively single interface
    λ = 1.0

    # With zero thickness, should get zero reflectance (matched media)
    Rs, Rp, Ts, Tp = airy(n0, nf, nf, 0.0, λ)
    Rs_expected, Rp_expected = fresnel(0.0, n0, nf)
    @test Rs ≈ Rs_expected atol=1e-10
    @test Rp ≈ Rp_expected atol=1e-10

    # Test 2: Energy conservation for lossless media
    n0 = 1.0
    nf = 1.5
    ns = 1.3
    d = 0.25
    λ = 1.0

    Rs, Rp, Ts, Tp = airy(n0, nf, ns, d, λ)
    @test Rs + Ts ≈ 1.0 atol=1e-10
    @test Rp + Tp ≈ 1.0 atol=1e-10
    @test Rs >= 0.0
    @test Rp >= 0.0
    @test Ts >= 0.0
    @test Tp >= 0.0

    # Test 3: Quarter-wave condition (n*d = λ/4)
    # For quarter-wave film: R = ((n0*ns - nf²)/(n0*ns + nf²))²
    n0 = 1.0
    nf = 1.38  # typical AR coating index
    ns = 1.52  # glass
    λ = 0.55
    d = λ / (4 * nf)  # quarter-wave thickness

    Rs, Rp, Ts, Tp = airy(n0, nf, ns, d, λ)

    # Analytical quarter-wave reflectance at normal incidence
    R_qw = ((n0 * ns - nf^2) / (n0 * ns + nf^2))^2
    @test Rs ≈ R_qw atol=1e-8
    @test Rp ≈ R_qw atol=1e-8

    # Test 4: Half-wave condition (n*d = λ/2)
    # Half-wave film is "absentee" - reflectance equals bare substrate
    d_half = λ / (2 * nf)
    Rs, Rp, Ts, Tp = airy(n0, nf, ns, d_half, λ)
    Rs_bare, Rp_bare = fresnel(0.0, n0, ns)
    @test Rs ≈ Rs_bare atol=1e-8
    @test Rp ≈ Rp_bare atol=1e-8

    # Test 5: Oblique incidence - still conserves energy
    θ = 0.5  # ~28.6 degrees
    Rs, Rp, Ts, Tp = airy(n0, nf, ns, d, λ; θ=θ)
    @test Rs + Ts ≈ 1.0 atol=1e-8
    @test Rp + Tp ≈ 1.0 atol=1e-8
    @test Rs != Rp  # s and p differ at oblique incidence

    # Test 6: Perfect AR coating (nf = √(n0*ns))
    # At quarter-wave thickness, R should be zero
    nf_ideal = √(n0 * ns)
    d_ideal = λ / (4 * nf_ideal)
    Rs, Rp, Ts, Tp = airy(n0, nf_ideal, ns, d_ideal, λ)
    @test Rs ≈ 0.0 atol=1e-10
    @test Rp ≈ 0.0 atol=1e-10
    @test Ts ≈ 1.0 atol=1e-10
    @test Tp ≈ 1.0 atol=1e-10
end

@testset "airy vs TMM" begin
    # Validate TMM against Airy formula for single thin film
    λ = 1.0
    λs = [λ, 1.1]
    n_air = 1.0
    n_film = 1.5
    n_sub = 1.3
    d = λ / (4 * n_film)

    # Build TMM layers
    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d)
    substrate = Layer(λs, fill(n_sub, length(λs)), zeros(length(λs)), 0.0)
    layers = [air, film, substrate]

    # TMM calculation
    Tpp_tmm, Tss_tmm, Rpp_tmm, Rss_tmm = transfer(λ, layers)

    # Airy calculation
    Rs_airy, Rp_airy, Ts_airy, Tp_airy = airy(n_air, n_film, n_sub, d, λ)

    # Compare - note TMM returns (Tpp, Tss, Rpp, Rss) order
    @test Rss_tmm ≈ Rs_airy atol=1e-6
    @test Rpp_tmm ≈ Rp_airy atol=1e-6
    @test Tss_tmm ≈ Ts_airy atol=1e-6
    @test Tpp_tmm ≈ Tp_airy atol=1e-6

    # Test at oblique incidence
    θ = 0.3
    Tpp_tmm, Tss_tmm, Rpp_tmm, Rss_tmm = transfer(λ, layers; θ=θ)
    Rs_airy, Rp_airy, Ts_airy, Tp_airy = airy(n_air, n_film, n_sub, d, λ; θ=θ)

    @test Rss_tmm ≈ Rs_airy atol=1e-6
    @test Rpp_tmm ≈ Rp_airy atol=1e-6
    @test Tss_tmm ≈ Ts_airy atol=1e-6
    @test Tpp_tmm ≈ Tp_airy atol=1e-6
end
