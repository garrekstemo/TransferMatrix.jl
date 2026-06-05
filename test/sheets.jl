using Test
using LinearAlgebra
using StaticArrays
using RefractiveIndex
using TransferMatrix

const Z0_test = sqrt(TransferMatrix.μ_0 / TransferMatrix.ε_0)

@testset "Sheet constructors" begin
    # Scalar constant -> diagonal isotropic tensor in SI Siemens
    s = TransferMatrix.Sheet(2.0e-4 + 1.0e-4im)
    σ = s.conductivity(1.0)
    @test σ isa SMatrix{2,2,ComplexF64}
    @test σ[1,1] == 2.0e-4 + 1.0e-4im
    @test σ[2,2] == 2.0e-4 + 1.0e-4im
    @test σ[1,2] == 0
    @test σ[2,1] == 0

    # Scalar function of λ
    sf = TransferMatrix.Sheet(λ -> 1.0e-4 / λ)
    @test sf.conductivity(2.0)[1,1] ≈ 5.0e-5
    @test sf.conductivity(2.0)[1,1] == sf.conductivity(2.0)[2,2]

    # Keyword anisotropic tensor (constants and functions mixed)
    st = TransferMatrix.Sheet(; xx = 1.0e-4, yy = 3.0e-4, xy = λ -> 1.0e-5)
    M = st.conductivity(1.0)
    @test M[1,1] ≈ 1.0e-4
    @test M[2,2] ≈ 3.0e-4
    @test M[1,2] ≈ 1.0e-5
    @test M[2,1] == 0

    # n,k + thickness -> σ = -i ω ε₀ d (n²-1), isotropic; agrees with closed form
    n = 4.0 + 0.2im
    d = 6.5e-4   # μm (~0.65 nm)
    λ = 0.6
    sheet_nk = TransferMatrix.Sheet(λ0 -> n, d)
    σ_expected = -im * (2π * TransferMatrix.c_0 / λ) * TransferMatrix.ε_0 * d * (n^2 - 1)
    @test sheet_nk.conductivity(λ)[1,1] ≈ σ_expected
    @test sheet_nk.conductivity(λ)[2,2] ≈ σ_expected
    @test sheet_nk.conductivity(λ)[1,2] == 0

    # In-plane anisotropic n,k
    sheet_aniso = TransferMatrix.Sheet(λ0 -> 4.0 + 0.0im, λ0 -> 5.0 + 0.0im, d)
    σx = -im * (2π * TransferMatrix.c_0 / λ) * TransferMatrix.ε_0 * d * ((4.0+0im)^2 - 1)
    σy = -im * (2π * TransferMatrix.c_0 / λ) * TransferMatrix.ε_0 * d * ((5.0+0im)^2 - 1)
    @test sheet_aniso.conductivity(λ)[1,1] ≈ σx
    @test sheet_aniso.conductivity(λ)[2,2] ≈ σy
end

@testset "sheet_matrix G structure" begin
    # G = G_phys^-1 in basis (Ex, Ey, H_y, -Hx): identity + (+σ̃) in rows 3,4.
    σxx = 1.0e-4 + 2.0e-4im
    σyy = 3.0e-4 - 1.0e-4im
    σxy = 0.5e-4im
    σyx = -0.2e-4im
    sheet = TransferMatrix.Sheet(; xx = σxx, yy = σyy, xy = σxy, yx = σyx)
    G = TransferMatrix.sheet_matrix(sheet, 1.0)

    @test G isa SMatrix{4,4,ComplexF64}
    # Top-left 2×2 identity, tangential E continuous
    @test G[1,1] == 1 && G[2,2] == 1
    @test G[1,2] == 0 && G[1,3] == 0 && G[1,4] == 0
    @test G[2,1] == 0 && G[2,3] == 0 && G[2,4] == 0
    # H rows carry +σ̃ = +Z₀ σ
    @test G[3,1] ≈ Z0_test * σxx
    @test G[3,2] ≈ Z0_test * σxy
    @test G[4,1] ≈ Z0_test * σyx
    @test G[4,2] ≈ Z0_test * σyy
    @test G[3,3] == 1 && G[4,4] == 1
    @test G[3,4] == 0 && G[4,3] == 0

    # Zero conductivity -> identity (no-op interface)
    G0 = TransferMatrix.sheet_matrix(TransferMatrix.Sheet(0.0 + 0.0im), 1.0)
    @test G0 ≈ SMatrix{4,4,ComplexF64}(I)
end
