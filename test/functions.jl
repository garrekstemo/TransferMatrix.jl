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


# @testset "poynting" begin

#     M = Array{ComplexF64}(undef, 6, 6)
#     ε = Diagonal([1, 1, 1])
#     μ = Diagonal([1, 1, 1])
#     M[1:3, 1:3] = ε
#     M[4:6, 4:6] = μ
    
#     a = TransferMatrix.construct_a(0., M)
#     Δ = TransferMatrix.construct_Δ(0., M, a)
#     q, Ψ = eigen(Δ)

#     # This just makes the elements of Ψ all 1 or -1
#     # for convenience when calculating by hand.
#     Ψ ./= √2 / 2 

#     # Ψ now looks like the following:
#     # Ψ = [1  0  1  0;
#     #     -1  0  1  0;
#     #      0  1  0  1;
#     #      0 -1  0  1]

#     # The relevant indices of the matrix, a, above are zero,
#     # so we give them arbitrary values here to check that
#     # all elements of the Poynting vector are calculated:
    
#     a[3,1] = 1
#     a[3,2] = 2
#     a[3,4] = 3
#     a[3,5] = 4
    
#     a[6,1] = 5
#     a[6,2] = 6
#     a[6,4] = 7
#     a[6,5] = 8
              
#     S = TransferMatrix.poynting(Ψ, a)

#     # Calculated by hand with the above a and Ψ following Passler et al. 2017.
#     S_test = ComplexF64[
#         -3  13 -5  -1;
#          3  5  -13  1;
#         -1 -1   1   1]

#     @test isapprox(S, S_test, atol=1e-13)
# end