using Test
using TransferMatrix

@testset "dispersion models" begin

    @testset "_to_eV seam is a no-op for Real" begin
        @test TransferMatrix._to_eV(2.0) === 2.0
        @test TransferMatrix._to_eV(9) === 9
    end

    @testset "photon energy conversion" begin
        # E[eV] = hc/λ; HC_EV_UM ≈ 1.2398 eV·μm
        @test TransferMatrix.HC_EV_UM ≈ 1.239841984 atol = 1e-6
        @test TransferMatrix._photon_energy_eV(TransferMatrix.HC_EV_UM / 2.0) ≈ 2.0
    end

    @testset "Drude: lossless below plasma is reflective" begin
        # ε_inf=1, ω_p=9 eV, γ=0 at E=2 eV → ε = 1 - 81/4 = -19.25, n purely imaginary
        n_drude = drude(9.0, 0.0)
        λ = TransferMatrix.HC_EV_UM / 2.0          # E = 2 eV
        n = n_drude(λ)
        @test real(n) ≈ 0.0 atol = 1e-9
        @test abs(imag(n)) ≈ sqrt(19.25)
        # purely imaginary index → unit reflectance at a single interface
        Rs, Rp = fresnel(0.0, 1.0, n)
        @test Rs ≈ 1.0 atol = 1e-12
        @test Rp ≈ 1.0 atol = 1e-12
    end

    @testset "Drude: damping makes ε″ > 0 (absorption, exp(-iωt))" begin
        n = drude(9.0, 0.3)(TransferMatrix.HC_EV_UM / 2.0)
        @test imag(n) > 0.0          # n″ > 0
        @test real(n) ≥ 0.0
    end

    @testset "Lorentz: static limit ε(0) = ε_inf + Δε" begin
        # ε_inf=2, single oscillator (ω_0=3, Δε=4, γ=0.1); as ω→0, ε→6, n→√6
        n_lor = lorentz(3.0, 4.0, 0.1; ε_inf=2.0)
        λ_long = TransferMatrix.HC_EV_UM / 1e-6     # E ≈ 1e-6 eV ≈ 0
        n = n_lor(λ_long)
        @test real(n) ≈ sqrt(6.0) atol = 1e-3
        @test imag(n) ≈ 0.0 atol = 1e-3
    end

    @testset "Lorentz: multi-oscillator static limit sums Δε" begin
        # ε_inf=1.5 + Δε of 1.0 and 0.5 → ε(0)=3.0, n→√3
        n_lor = lorentz([(2.0, 1.0, 0.1), (4.0, 0.5, 0.2)]; ε_inf=1.5)
        n = n_lor(TransferMatrix.HC_EV_UM / 1e-6)
        @test real(n) ≈ sqrt(3.0) atol = 1e-3
        @test imag(n) ≈ 0.0 atol = 1e-3
    end

    @testset "Lorentz: lossless, far below resonance is real" begin
        # γ=0, evaluate well below ω_0 → ε real > 0 → n real
        n = lorentz(5.0, 1.0, 0.0; ε_inf=2.0)(TransferMatrix.HC_EV_UM / 1.0)  # E=1 eV ≪ 5 eV
        @test imag(n) ≈ 0.0 atol = 1e-12
        @test real(n) > sqrt(2.0)        # resonance below adds to ε∞
    end

    @testset "Lorentz: single-arg delegates to one-oscillator vector" begin
        λ = TransferMatrix.HC_EV_UM / 2.0
        @test lorentz(3.0, 4.0, 0.1; ε_inf=2.0)(λ) ≈
              lorentz([(3.0, 4.0, 0.1)]; ε_inf=2.0)(λ)
    end

end
