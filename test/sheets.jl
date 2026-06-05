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

@testset "sheet σ=0 regression (transfer)" begin
    n1, n2, nf = 1.0, 1.5, 2.0
    d = 0.1
    air = Layer(λ -> complex(n1), 0.0)
    film = Layer(λ -> complex(nf), d)
    sub = Layer(λ -> complex(n2), 0.0)
    layers = [air, film, sub]

    base = transfer(0.6, layers; θ = 0.2)
    zero_sheet = Dict(1 => TransferMatrix.Sheet(0.0 + 0.0im))
    withσ0 = transfer(0.6, layers; θ = 0.2, sheets = zero_sheet)

    @test isapprox(base.Rpp, withσ0.Rpp; atol = 1e-12)
    @test isapprox(base.Rss, withσ0.Rss; atol = 1e-12)
    @test isapprox(base.Tpp, withσ0.Tpp; atol = 1e-12)
    @test isapprox(base.Tss, withσ0.Tss; atol = 1e-12)
end

@testset "analytic conductive interface (pins G sign)" begin
    # Single sheet between media n1 | n2; compare R to closed-form conductive Fresnel.
    n1, n2 = 1.0, 1.5
    σ_s = 2.0e-4 + 1.0e-4im
    g = Z0_test * σ_s                      # dimensionless σ̃
    air = Layer(λ -> complex(n1), 0.0)
    sub = Layer(λ -> complex(n2), 0.0)
    layers = [air, sub]
    sheets = Dict(1 => TransferMatrix.Sheet(σ_s))

    for θ in (0.0, π/6)
        c1 = cos(θ)
        c2 = sqrt(1 - (n1 / n2 * sin(θ))^2)        # real (n1 < n2)
        rs = (n1*c1 - n2*c2 - g) / (n1*c1 + n2*c2 + g)
        rp = (n2*c1 - n1*c2 + g*c1*c2) / (n2*c1 + n1*c2 + g*c1*c2)
        res = transfer(1.0, layers; θ = θ, sheets = sheets)
        @test isapprox(res.Rss, abs2(rs); atol = 1e-8)
        @test isapprox(res.Rpp, abs2(rp); atol = 1e-8)
    end
end
