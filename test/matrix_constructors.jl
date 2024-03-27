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

@testset "dynamical_matrix" begin
    γ = ComplexF64[
         11 12 13;
         21 22 23;
         31 32 33;
         41 42 43]
    q = ComplexF64[1, 2, 3, 4]
    ξ = 1.0 + 1.0im
    μ = 2.0

    A_test = ComplexF64[
        11 21 31 41
        12 22 32 42
     -2 - 13im  19 - 23im    60 - 33im  121 - 43im
        12 44 96 168
    ]
    A_test[3, :] ./= μ
    A_test[4, :] ./= μ

    A = TransferMatrix.dynamical_matrix(ξ, q, γ, μ)

    for (i, col) in enumerate(eachrow(A))
        @test A[:, i] == A_test[:, i]
    end
end

@testset "propagation_matrix" begin
    q = [1., 2., 3., 4.]
    d = π / 2
    ω = c_0
    P_test = Diagonal([-1.0im, -1.0, 1.0im, 1.0])
    P = TransferMatrix.propagation_matrix(ω, q)

    @test P(d) ≈ P_test
end


@testset "calculate_γ" begin

    μ = 1.0

    ξ1 = 0.0
    ξ2 = 1.0
    ξ3 = 1.0 + 1im

    ε1 = I
    ε2 = fill(1, (3, 3))
    ε3 = [3 1 1 ; 1 2 1 ; 1 1 2]
    ε4 = [1 1 1 ; 1 1 1 ; 1 1 2]
    ε5 = [1 + 1im   1        1;
            1     1 + 1im    1;
            1       1      2 + 2im]

    q1 = [1., 1., 1., 1.]
    q2 = [1., 2., 1., 2.]
    q3 = [1. + 1im, 1. + 1im, 1. + 1im, 1. + 1im]
    q4 = [1. + 1im, 2. + 1im, 1. + 1im, 2. + 1im]

    γ1 = TransferMatrix.calculate_γ(ξ1, q1, ε1, μ)
    γ2 = TransferMatrix.calculate_γ(ξ1, q1, ε2, μ)
    γ3 = TransferMatrix.calculate_γ(ξ2, q2, ε3, μ)
    γ4 = TransferMatrix.calculate_γ(ξ2, q1, ε4, μ)
    γ5 = TransferMatrix.calculate_γ(ξ3, q3, ε5, μ)
    γ6 = TransferMatrix.calculate_γ(ξ3, q4, ε5, μ)

    γ1_test = ComplexF64[
        1 0 0;
        0 1 0;
       -1 0 0;
        0 1 0
        ]

    γ2_test = ComplexF64[
        1  0 -1;
        0  1 -1;
       -1  0  1;
        0  1 -1
        ]

    γ3_test = ComplexF64[
        1.0 -1.0 -1.0;
       -0.2  1.0 -0.4;
       -1.0  1.0  3.0;
       -0.2  1.0 -0.4
    ]

    γ4_test = ComplexF64[
        1.0  0.0 -2.0;
        0.0  1.0 -1.0;
       -1.0  0.0  2.0;
        0.0  1.0 -1.0
    ]

    γ5_test = ComplexF64[
        1.0  0.0 -0.5 - 1im;
        0.0  1.0 -0.5;
       -1.0  0.0  0.5 + 1im;
        0.0  1.0 -0.5
    ]

    γ6_test = ComplexF64[       
          1.0+0.0im         -0.351351-0.108108im  -0.324324-0.945946im
    -0.166154+0.00923077im        1.0+0.0im           -0.32+0.24im
         -1.0+0.0im          0.351351+0.108108im   0.675676+1.05405im
    -0.166154+0.00923077im        1.0+0.0im           -0.32+0.24im
    ]

    for m in 1:4
        γ2_test[m, :] /= √(γ2_test[m, :] ⋅ γ2_test[m, :]')
        γ3_test[m, :] /= √(γ3_test[m, :] ⋅ γ3_test[m, :]')
        γ4_test[m, :] /= √(γ4_test[m, :] ⋅ γ4_test[m, :]')
        γ5_test[m, :] /= √(γ5_test[m, :] ⋅ γ5_test[m, :]')
        γ6_test[m, :] /= √(γ6_test[m, :] ⋅ γ6_test[m, :]')
    end

    @test γ1 == γ1_test
    @test γ2 == γ2_test
    @test isapprox(γ3, γ3_test)
    @test isapprox(γ4, γ4_test)
    @test isapprox(γ5, γ5_test)
    @test isapprox(γ6, γ6_test, atol=1e-5)
end