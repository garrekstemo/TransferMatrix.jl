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

end
