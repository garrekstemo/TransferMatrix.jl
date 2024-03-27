using Test
using LinearAlgebra
using StaticArrays
using RefractiveIndex
using TransferMatrix

const c_0 = 299792458

@testset "construct_a" begin

    # Test orthorhombic crystal with principal axes
    # parallel to x, y, z. M is constant and diagonal.
    # The only nonzero coefficients are
    
    # a[3,5] = -ξ / M[3, 3]
    # a[6,2] =  ξ / M[6, 6]
    
    # Berreman, Optics in Stratefied and Anisotropic Media, 1972
    # DOI: 10.1364/JOSA.62.000502

    ξ = 15.0 + 15im

    ε = Diagonal([1 + 1im, 1 + 1im, 1 + 1im])
    μ = Diagonal([1 + 0im, 1 + 0im, 1 + 0im])
    ε[3,3] = 3.0 + 0im
    μ[3,3] = 0.0 + 5im

    M = @MMatrix zeros(ComplexF64, 6, 6)
    M[1:3, 1:3] = ε
    M[4:6, 4:6] = μ
    
    a = TransferMatrix.construct_a(ξ, M)

    to_subtract = @MMatrix zeros(ComplexF64, 6, 6)
    to_subtract[3,5] = -5.0 - 5im
    to_subtract[6,2] = 3.0 - 3im

    b = a - to_subtract
    test_against = zeros(ComplexF64, 6, 6)

    @test isapprox(b, zeros(ComplexF64, 6, 6), atol=1e-15)
end

@testset "construct_Δ" begin
 
    # Continue the test from `construct_a`.
    # For an orthorhombic crystal described above, the only nonzero elements of Δ are
    # Δ[2,1] = M[1,1] = ε[1,1]
    # Δ[4,3] = M[2,2] - ξ^2 / M[6,6] = ε[2,2] - ξ^2 / μ[3,3]
    # Δ[3,4] = M[4,4] = μ[1,1]
    # Δ[1,2] = M[5,5] - ξ^2 / M[3,3] = μ[2,2] - ξ^2 / ε[3,3]

    ξ = 0.

    ε = Diagonal([1 + 1im, 1 + 1im, 1 + 1im])
    μ = Diagonal([1 + 0im, 1 + 0im, 1 + 0im])
    ε[3,3] = 3.0 + 0im
    μ[3,3] = 0.0 + 5im

    M = zeros(ComplexF64, 6, 6)
    M[1:3, 1:3] = ε
    M[4:6, 4:6] = μ
    
    a = TransferMatrix.construct_a(ξ, M)
    Δ = TransferMatrix.construct_Δ(ξ, M, a)

    Δ21 = ε[1,1]
    Δ43 = ε[2,2] - ξ^2 / μ[3,3]
    Δ34 = μ[1,1]
    Δ12 = μ[2,2] - ξ^2 / ε[3,3]
    Δ_squared = Diagonal([Δ12 * Δ21, Δ12 * Δ21, Δ34 * Δ43, Δ34 * Δ43])
    Δ_cubed = Δ^3

    @test Δ^2 == Δ_squared
    @test Δ_cubed[1,1] == 0.0 + 0im
    @test Δ_cubed[1,2] == Δ12^2 * Δ21
end

@testset "construct_M" begin
    ε_i = (1.0 + 2.0im)^2
    μ_i = 1.0 + 0.0im
    ε = Diagonal([ε_i, ε_i, ε_i])
    μ = Diagonal(fill(μ_i, 3))
    ρ1 = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    ρ2 = [9.0 8.0 7.0; 6.0 5.0 4.0; 3.0 2.0 1.0]

    M1 = TransferMatrix.construct_M(ε, μ, ρ1, ρ2)
    M1_true = [
        ε_i 0.0 0.0 1.0 2.0 3.0;
        0.0 ε_i 0.0 4.0 5.0 6.0;
        0.0 0.0 ε_i 7.0 8.0 9.0;
        9.0 8.0 7.0 μ_i 0.0 0.0;
        6.0 5.0 4.0 0.0 μ_i 0.0;
        3.0 2.0 1.0 0.0 0.0 μ_i
    ]

    M2 = TransferMatrix.construct_M(ε, μ)
    M2_true = [
        ε_i 0.0 0.0 0.0 0.0 0.0;
        0.0 ε_i 0.0 0.0 0.0 0.0;
        0.0 0.0 ε_i 0.0 0.0 0.0;
        0.0 0.0 0.0 μ_i 0.0 0.0;
        0.0 0.0 0.0 0.0 μ_i 0.0;
        0.0 0.0 0.0 0.0 0.0 μ_i
    ]

    M3 = TransferMatrix.construct_M(ε)
    M3_true = [
        ε_i 0.0 0.0 0.0 0.0 0.0;
        0.0 ε_i 0.0 0.0 0.0 0.0;
        0.0 0.0 ε_i 0.0 0.0 0.0;
        0.0 0.0 0.0 1.0 0.0 0.0;
        0.0 0.0 0.0 0.0 1.0 0.0;
        0.0 0.0 0.0 0.0 0.0 1.0
    ]

    @test M1 == M1_true
    @test M2 == M2_true
    @test M3 == M3_true
end