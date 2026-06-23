# Tests for circular-polarization basis output (issue #84).

@testset "circular basis-change helper" begin
    # Helicity flip on reflection: r_lin = diag(0.2, -0.2) -> purely off-diagonal.
    r_lin = SMatrix{2,2,ComplexF64}(0.2, 0.0, 0.0, -0.2)  # [0.2 0; 0 -0.2]
    rc = TransferMatrix._jones_to_circular(r_lin)
    @test isapprox(abs(rc[1, 1]), 0.0; atol=1e-12)   # LL
    @test isapprox(abs(rc[2, 2]), 0.0; atol=1e-12)   # RR
    @test isapprox(abs(rc[1, 2]), 0.2; atol=1e-12)   # LR
    @test isapprox(abs(rc[2, 1]), 0.2; atol=1e-12)   # RL

    # No flip on transmission: t_lin = diag(0.8, 0.8) -> purely diagonal.
    t_lin = SMatrix{2,2,ComplexF64}(0.8, 0.0, 0.0, 0.8)
    tc = TransferMatrix._jones_to_circular(t_lin)
    @test isapprox(abs(tc[1, 2]), 0.0; atol=1e-12)
    @test isapprox(abs(tc[2, 1]), 0.0; atol=1e-12)
    @test isapprox(tc[1, 1], 0.8 + 0im; atol=1e-12)
    @test isapprox(tc[2, 2], 0.8 + 0im; atol=1e-12)

    # Over-correction guard: ideal mirror r_lin = -I stays DIAGONAL (-I),
    # because -I is invariant under any unitary similarity.
    mirror = SMatrix{2,2,ComplexF64}(-1.0, 0.0, 0.0, -1.0)
    mc = TransferMatrix._jones_to_circular(mirror)
    @test isapprox(abs(mc[1, 2]), 0.0; atol=1e-12)
    @test isapprox(abs(mc[2, 1]), 0.0; atol=1e-12)
    @test isapprox(mc[1, 1], -1.0 + 0im; atol=1e-12)

    # Closed form for diagonal (isotropic) input diag(a,b):
    # [ (a+b)/2  (a-b)/2; (a-b)/2  (a+b)/2 ].
    a, b = 0.3 + 0im, -0.5 + 0im
    diagj = SMatrix{2,2,ComplexF64}(a, 0.0, 0.0, b)
    jc = TransferMatrix._jones_to_circular(diagj)
    @test isapprox(jc[1, 1], (a + b) / 2; atol=1e-12)
    @test isapprox(jc[2, 2], (a + b) / 2; atol=1e-12)
    @test isapprox(jc[1, 2], (a - b) / 2; atol=1e-12)
    @test isapprox(jc[2, 1], (a - b) / 2; atol=1e-12)

    # Unitarity: total |J|^2 is basis-invariant.
    J = SMatrix{2,2,ComplexF64}(0.2 + 0.1im, -0.05im, 0.3, 0.4 - 0.2im)
    @test isapprox(sum(abs2, TransferMatrix._jones_to_circular(J)), sum(abs2, J); atol=1e-12)

    # Round-trip: C * (C^{-1} J C) * C^{-1} recovers J.
    back = TransferMatrix._C_CIRC * TransferMatrix._jones_to_circular(J) * TransferMatrix._C_CIRC_INV
    @test isapprox(back, J; atol=1e-12)
end

@testset "circular result assembly" begin
    # Synthetic isotropic single-interface (air->glass) amplitudes:
    # r = (rpp, rps, rss, rsp), t = (tpp, tps, tsp, tss).
    r = SVector{4,ComplexF64}(0.2, 0.0, -0.2, 0.0)
    t = SVector{4,ComplexF64}(0.8, 0.0, 0.0, 0.8)
    Tpp = 0.96  # Poynting transmittance (N = (0.96+0.96)/(0.64+0.64) = 1.5)
    Tss = 0.96

    res = TransferMatrix._circular_result(r, t, Tpp, Tss)
    @test res isa TransferMatrix.CircularTransferResult{Float64}

    # Reflection flips: diagonal ~0, cross = 0.04.
    @test isapprox(res.Rrr, 0.0; atol=1e-12)
    @test isapprox(res.Rll, 0.0; atol=1e-12)
    @test isapprox(res.Rrl, 0.04; atol=1e-12)
    @test isapprox(res.Rlr, 0.04; atol=1e-12)

    # Transmission preserves (diagonal), scaled by N = 1.5: 1.5 * 0.64 = 0.96.
    @test isapprox(res.Trr, 0.96; atol=1e-12)
    @test isapprox(res.Tll, 0.96; atol=1e-12)
    @test isapprox(res.Trl, 0.0; atol=1e-12)
    @test isapprox(res.Tlr, 0.0; atol=1e-12)

    # N guard: zero transmission amplitude -> all circular T are 0 (not NaN).
    res0 = TransferMatrix._circular_result(r, SVector{4,ComplexF64}(0, 0, 0, 0), 0.0, 0.0)
    @test res0.Trr == 0.0
    @test res0.Tll == 0.0
    @test !isnan(res0.Trr)
end
