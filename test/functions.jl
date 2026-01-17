# Test other functions in TransferMatrix.jl

@testset "abs_ratio" begin
    S = ComplexF64[
        1  1  0  0;
        1  2  0  0;
        0  0  0  0
        ]

    a = TransferMatrix.abs_ratio(S[1, 1], S[2, 1])
    b = TransferMatrix.abs_ratio(S[1, 2], S[2, 2])
    c = TransferMatrix.abs_ratio(S[1, 3], S[2, 3])
    d = TransferMatrix.abs_ratio(S[1, 4], S[2, 4])

    @test isapprox(a, b) == false
    @test a > b
    @test isapprox(c, d) == false
    @test isnan(c) == true
    @test isnan(d) == true
end


@testset "evaluate_birefringence" begin
    Ψ1 = ComplexF64[
        1 2 1 2;
        0 0 0 0;
        2 1 2 1;
        0 0 0 0
        ]

    # Evaluate 2/5 > 1/5 in reflected modes
    S = ComplexF64[
        1 1 1 2;
        1 1 2 1;
        0 0 0 0
        ]

    # Setting elements of the Poynting vector to zero
    # forces the electric field to be evaluated.
    S1 = zeros(Int64, 3, 4)

    t1, r1 = TransferMatrix.evaluate_birefringence(Ψ1, S, [1, 2], [3, 4])
    t2, r2 = TransferMatrix.evaluate_birefringence(Ψ1, S, [2, 1], [3, 4])
    t3, r3 = TransferMatrix.evaluate_birefringence(Ψ1, S1, [1, 2], [3, 4])
    t4, r4 = TransferMatrix.evaluate_birefringence(Ψ1, S1, [2, 1], [4, 3])

    # transmitted modes not reversed in either case
    @test t1 == [1, 2]
    @test t2 == [2, 1]

    # reflected modes reversed
    @test r1 == [4, 3]
    @test r1 != [3, 4]

    # both transmitted and reflected modes reversed
    @test t3 == [2, 1]
    @test r3 == [4, 3]

    # neither are reversed
    @test t4 == [2, 1]
    @test r4 == [4, 3]
end


@testset "poynting(Ψ, a)" begin

    M = Array{ComplexF64}(undef, 6, 6)
    ε = Diagonal([1, 1, 1])
    μ = Diagonal([1, 1, 1])
    M[1:3, 1:3] = ε
    M[4:6, 4:6] = μ
    
    a = Matrix(TransferMatrix.construct_a(0., M))

    # Use a fixed Ψ to keep the expected Poynting result stable.
    Ψ = ComplexF64[
        1  0  1  0;
       -1  0  1  0;
        0  1  0  1;
        0 -1  0  1
    ]

    # The relevant indices of the matrix, a, above are zero,
    # so we give them arbitrary values here to check that
    # all elements of the Poynting vector are calculated:
    
    a[3,1] = 1
    a[3,2] = 2
    a[3,4] = 3
    a[3,5] = 4
    
    a[6,1] = 5
    a[6,2] = 6
    a[6,4] = 7
    a[6,5] = 8
              
    S = TransferMatrix.poynting(Ψ, a)

    # Calculated by hand with the above a and Ψ following Passler et al. 2017.
    S_test = ComplexF64[
        -3  13 -5  -1;
         3  5  -13  1;
        -1 -1   1   1]

    @test isapprox(S, S_test, atol=1e-13)
end

@testset "poynting(ξ, q_in, q_out, γ_in, γ_out, t_coefs, r_coefs)" begin
    ξ = 0.0
    q_in = SVector(1.0, 1.0, 1.0, 1.0)
    q_out = SVector(1.0, 1.0, 1.0, 1.0)
    γ_in = @SMatrix [1.0 0.0 0.0;
                     0.0 1.0 0.0;
                     0.0 0.0 1.0;
                     1.0 0.0 0.0]
    γ_out = γ_in
    t_coefs = zeros(4)
    r_coefs = zeros(4)

    S = TransferMatrix.poynting(ξ, q_in, q_out, γ_in, γ_out, t_coefs, r_coefs)
    @test norm(S.in_p) > 0
    @test norm(S.in_s) > 0
    @test isapprox(S.out_p, SVector(0.0, 0.0, 0.0); atol=1e-12)
    @test isapprox(S.out_s, SVector(0.0, 0.0, 0.0); atol=1e-12)
    @test isapprox(S.refl_p, SVector(0.0, 0.0, 0.0); atol=1e-12)
    @test isapprox(S.refl_s, SVector(0.0, 0.0, 0.0); atol=1e-12)
end

@testset "calculate_q sorting" begin
    Δ = Diagonal([1.0, 2.0, -1.0, -2.0])
    a = zeros(ComplexF64, 6, 6)
    q, S = TransferMatrix.calculate_q(Matrix(Δ), a)
    @test length(q) == 4
    @test all(q[1:2] .> 0)
    @test all(q[3:4] .< 0)
end

@testset "calculate_q complex sorting" begin
    Δ = Diagonal(ComplexF64[1 + 1im, 2 + 2im, 3 - 1im, 4 - 2im])
    a = zeros(ComplexF64, 6, 6)
    q, S = TransferMatrix.calculate_q(Matrix(Δ), a)
    @test length(q) == 4
    @test all(imag.(q[1:2]) .>= 0)
    @test all(imag.(q[3:4]) .<= 0)
end

@testset "calculate_q mode sorting error" begin
    a = zeros(ComplexF64, 6, 6)

    # All real eigenvalues positive → 4 transmitted, 0 reflected
    Δ_all_positive = Diagonal([1.0, 2.0, 3.0, 4.0])
    @test_throws ArgumentError TransferMatrix.calculate_q(Matrix(Δ_all_positive), a)

    # All real eigenvalues negative → 0 transmitted, 4 reflected
    Δ_all_negative = Diagonal([-1.0, -2.0, -3.0, -4.0])
    @test_throws ArgumentError TransferMatrix.calculate_q(Matrix(Δ_all_negative), a)

    # 3 positive, 1 negative → 3 transmitted, 1 reflected
    Δ_three_one = Diagonal([1.0, 2.0, 3.0, -1.0])
    @test_throws ArgumentError TransferMatrix.calculate_q(Matrix(Δ_three_one), a)

    # Complex: all positive imaginary → 4 transmitted, 0 reflected
    Δ_complex_all_pos = Diagonal(ComplexF64[1 + 1im, 2 + 2im, 3 + 3im, 4 + 4im])
    @test_throws ArgumentError TransferMatrix.calculate_q(Matrix(Δ_complex_all_pos), a)

    # Complex: all negative imaginary → 0 transmitted, 4 reflected
    Δ_complex_all_neg = Diagonal(ComplexF64[1 - 1im, 2 - 2im, 3 - 3im, 4 - 4im])
    @test_throws ArgumentError TransferMatrix.calculate_q(Matrix(Δ_complex_all_neg), a)

    # Verify error message contains eigenvalue info
    try
        TransferMatrix.calculate_q(Matrix(Δ_all_positive), a)
        @test false  # Should not reach here
    catch e
        @test e isa ArgumentError
        @test occursin("Mode sorting failed", e.msg)
        @test occursin("4 transmitted", e.msg)
        @test occursin("0 reflected", e.msg)
    end
end

@testset "dielectric_tensor" begin
    ε = TransferMatrix.dielectric_tensor(2.25, 2.0, 1.5)
    @test ε isa Diagonal
    @test ε[1,1] ≈ 2.25
    @test ε[2,2] ≈ 2.0
    @test ε[3,3] ≈ 1.5
    @test ε[1,2] == 0
    @test ε[1,3] == 0
    @test ε[2,3] == 0
end

@testset "permeability_tensor" begin
    μ = TransferMatrix.permeability_tensor(1.0, 2.0, 3.0)
    @test μ isa Diagonal
    @test μ[1,1] ≈ 1.0
    @test μ[2,2] ≈ 2.0
    @test μ[3,3] ≈ 3.0
    @test μ[1,2] == 0
end

@testset "construct_M overloads" begin
    # Test general construct_M with magnetoelectric coupling
    ε = Diagonal([2.25, 2.25, 2.25])
    μ = Diagonal([1.0, 1.0, 1.0])
    ρ1 = [0.1 0 0; 0 0.1 0; 0 0 0.1]
    ρ2 = [0.2 0 0; 0 0.2 0; 0 0 0.2]

    M = TransferMatrix.construct_M(ε, μ, ρ1, ρ2)
    @test size(M) == (6, 6)
    @test M[1:3, 1:3] == ε
    @test M[4:6, 4:6] == μ
    @test M[1:3, 4:6] == ρ1
    @test M[4:6, 1:3] == ρ2

    # Test specialized Diagonal method (isotropic, no magnetoelectric coupling)
    ε_diag = Diagonal(SVector{3, ComplexF64}(2.25, 2.25, 2.25))
    μ_diag = Diagonal(SVector{3, ComplexF64}(1.0, 1.0, 1.0))
    M_iso = TransferMatrix.construct_M(ε_diag, μ_diag)
    @test M_iso isa SMatrix{6,6}
    @test M_iso[1,1] ≈ 2.25
    @test M_iso[4,4] ≈ 1.0
    @test M_iso[1,4] == 0  # No magnetoelectric coupling
    @test M_iso[4,1] == 0

    # Test specialized SMatrix{3,3} method (rotated anisotropic)
    ε_full = @SMatrix [2.25+0im 0.1+0im 0.0+0im;
                       0.1+0im 2.0+0im 0.0+0im;
                       0.0+0im 0.0+0im 1.5+0im]
    M_aniso = TransferMatrix.construct_M(ε_full, μ_diag)
    @test M_aniso isa SMatrix{6,6}
    @test M_aniso[1,1] ≈ 2.25
    @test M_aniso[1,2] ≈ 0.1  # Off-diagonal element preserved
    @test M_aniso[2,1] ≈ 0.1
    @test M_aniso[4,4] ≈ 1.0
    @test M_aniso[1,4] == 0  # No magnetoelectric coupling
end

@testset "construct_a" begin
    # Test with isotropic material at normal incidence
    ε = Diagonal(SVector{3, ComplexF64}(2.25, 2.25, 2.25))
    μ = Diagonal(SVector{3, ComplexF64}(1.0, 1.0, 1.0))
    M = TransferMatrix.construct_M(ε, μ)

    ξ = 0.0  # Normal incidence
    a = TransferMatrix.construct_a(ξ, M)
    @test a isa SMatrix{6,6}

    # At normal incidence with isotropic material, most elements should be simple
    @test a[1,1] == 0  # First two rows are zero
    @test a[2,2] == 0

    # Test with oblique incidence
    ξ = 0.5
    a_oblique = TransferMatrix.construct_a(ξ, M)
    @test a_oblique isa SMatrix{6,6}
    # The a matrix has non-trivial values in rows 3 and 6
    @test !all(iszero, a_oblique[3,:])
end

@testset "construct_Δ" begin
    ε = Diagonal(SVector{3, ComplexF64}(2.25, 2.25, 2.25))
    μ = Diagonal(SVector{3, ComplexF64}(1.0, 1.0, 1.0))
    M = TransferMatrix.construct_M(ε, μ)

    ξ = 0.0
    a = TransferMatrix.construct_a(ξ, M)
    Δ = TransferMatrix.construct_Δ(ξ, M, a)

    @test Δ isa SMatrix{4,4}
    # For isotropic material at normal incidence, Δ should have specific structure
    # The eigenvalues should be ±n (refractive index)
    n = sqrt(2.25)
    eigvals = eigvals!(Matrix(Δ))
    @test any(v -> isapprox(abs(v), n, atol=1e-10), eigvals)
end

@testset "calculate_γ edge cases" begin
    # Test with degenerate q values (q[1] ≈ q[2])
    # This triggers the isapprox(q[1], q[2]) branch
    q_degen = SVector(1.5+0im, 1.5+0im, -1.5+0im, -1.5+0im)
    ε = Diagonal(SVector{3, ComplexF64}(2.25, 2.25, 2.25))
    ξ = 0.0
    μ = 1.0

    γ = TransferMatrix.calculate_γ(ξ, q_degen, ε, μ)
    @test γ isa SMatrix{4,3}
    # With degenerate eigenvalues, specific elements should be zero
    @test γ[1,2] ≈ 0
    @test γ[2,1] ≈ 0
    @test γ[3,2] ≈ 0
    @test γ[4,1] ≈ 0

    # Test singular denom_33 case (μ * ε[3,3] - ξ² ≈ 0)
    # This happens when ξ² = μ * ε[3,3]
    ε_singular = Diagonal(SVector{3, ComplexF64}(2.25, 2.25, 0.25))  # ε[3,3] = 0.25
    ξ_singular = 0.5  # ξ² = 0.25 = μ * ε[3,3] when μ = 1
    q_singular = SVector(1.0+0im, 0.5+0im, -1.0+0im, -0.5+0im)

    γ_singular = TransferMatrix.calculate_γ(ξ_singular, q_singular, ε_singular, μ)
    @test γ_singular isa SMatrix{4,3}
    # γ[:,3] elements should be zero due to singular_33 branch
    @test γ_singular[1,3] ≈ 0
    @test γ_singular[2,3] ≈ 0
    @test γ_singular[3,3] ≈ 0
    @test γ_singular[4,3] ≈ 0
end

@testset "dynamical_matrix" begin
    ξ = 0.0
    q = SVector(1.5+0im, 1.5+0im, -1.5+0im, -1.5+0im)
    γ = @SMatrix [1.0+0im 0.0+0im 0.0+0im;
                  0.0+0im 1.0+0im 0.0+0im;
                  -1.0+0im 0.0+0im 0.0+0im;
                  0.0+0im 1.0+0im 0.0+0im]
    μ = 1.0

    D = TransferMatrix.dynamical_matrix(ξ, q, γ, μ)
    @test D isa SMatrix{4,4}
    # First two rows should match γ transposed for columns 1 and 2
    @test D[1,1] ≈ γ[1,1]
    @test D[2,1] ≈ γ[1,2]
end

@testset "propagation_matrix" begin
    ω = 2π * c_0 / 500e-9  # 500 nm wavelength
    q = SVector(1.5+0im, 1.5+0im, -1.5+0im, -1.5+0im)

    P = TransferMatrix.propagation_matrix(ω, q)
    @test P isa TransferMatrix.PropagationMatrix

    # Test the callable interface
    z = 100e-9  # 100 nm thickness
    P_z = P(z)
    @test P_z isa Diagonal
    @test size(P_z) == (4, 4)

    # Propagation matrix elements should be complex exponentials with magnitude 1
    # for real q values
    for i in 1:4
        if imag(q[i]) == 0
            @test abs(P_z[i,i]) ≈ 1.0 atol=1e-10
        end
    end

    # Test that z=0 gives identity
    P_0 = P(0.0)
    @test P_0 ≈ Diagonal(ones(ComplexF64, 4))

    # Test that propagation is phase accumulation
    # P(z) = diag(exp(-iωqz/c))
    expected = exp.(-im * ω * q * z / c_0)
    for i in 1:4
        @test P_z[i,i] ≈ expected[i]
    end
end
