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
