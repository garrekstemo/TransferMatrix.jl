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
    # Poynting flux factors: unit-amplitude substrate-mode flux 1.5, unit
    # incident flux 1.0, so each |t|² converts with factor 1.5 (0.64 -> 0.96).
    f_out1 = f_out2 = 1.5
    f_in_p = f_in_s = 1.0

    res = TransferMatrix._circular_result(r, t, f_out1, f_out2, f_in_p, f_in_s)
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

    # Null guard: zero transmission amplitude -> all circular T are 0 (not NaN).
    res0 = TransferMatrix._circular_result(r, SVector{4,ComplexF64}(0, 0, 0, 0),
                                           f_out1, f_out2, f_in_p, f_in_s)
    @test res0.Trr == 0.0
    @test res0.Tll == 0.0
    @test !isnan(res0.Trr)
end

@testset "transfer basis=:circular single interface" begin
    # Air -> glass single interface at normal incidence (non-vacuum isotropic
    # substrate, so N = n2/n1 = 1.5).
    λs = [1.0, 2.0]
    air = Layer(λs, fill(1.0, length(λs)), zeros(length(λs)), 0.0)
    glass = Layer(λs, fill(1.5, length(λs)), zeros(length(λs)), 0.0)
    layers = [air, glass]
    λ = 1.5

    res = transfer(λ, layers; basis=:circular)
    @test res isa CircularTransferResult

    # Reflection flips helicity: diagonal ~0, cross = R_fresnel = 0.04.
    @test isapprox(res.Rrr, 0.0; atol=1e-10)
    @test isapprox(res.Rll, 0.0; atol=1e-10)
    @test isapprox(res.Rrl, 0.04; atol=1e-8)
    @test isapprox(res.Rlr, 0.04; atol=1e-8)
    @test isapprox(res.Rrl, res.Rlr; atol=1e-12)

    # Transmission preserves helicity (cross ~0), energy-conserving with N=1.5.
    @test isapprox(res.Trl, 0.0; atol=1e-10)
    @test isapprox(res.Tlr, 0.0; atol=1e-10)
    @test isapprox(res.Trr, 0.96; atol=1e-8)

    # Energy conservation per input polarization.
    @test isapprox(res.Rrr + res.Rlr + res.Trr + res.Tlr, 1.0; atol=1e-6)  # input R
    @test isapprox(res.Rll + res.Rrl + res.Tll + res.Trl, 1.0; atol=1e-6)  # input L
end

@testset "transfer linear path unchanged" begin
    λs = [1.0, 2.0]
    air = Layer(λs, fill(1.0, length(λs)), zeros(length(λs)), 0.0)
    glass = Layer(λs, fill(1.5, length(λs)), zeros(length(λs)), 0.0)
    layers = [air, glass]
    λ = 1.5

    default = transfer(λ, layers)
    explicit = transfer(λ, layers; basis=:linear)
    @test default isa TransferResult
    @test explicit isa TransferResult
    for f in (:Tpp, :Tss, :Tps, :Tsp, :Rpp, :Rss, :Rps, :Rsp)
        @test getfield(default, f) === getfield(explicit, f)
    end

    @test_throws ArgumentError transfer(λ, layers; basis=:bogus)
end

@testset "circular oblique isotropic symmetry and energy" begin
    λs = [1.0, 1.2]
    air = Layer(λs, fill(1.0, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(1.6, length(λs)), zeros(length(λs)), 0.25)
    layers = [air, film, air]  # vacuum substrate (N = 1)
    λ = 1.0
    θ = 0.3

    res = transfer(λ, layers; θ=θ, basis=:circular)
    # Isotropic media: Rrl=Rlr and Rrr=Rll even where Rpp != Rss.
    @test isapprox(res.Rrl, res.Rlr; atol=1e-10)
    @test isapprox(res.Rrr, res.Rll; atol=1e-10)
    # Vacuum substrate + lossless -> energy conserved per input polarization.
    @test isapprox(res.Rrr + res.Rlr + res.Trr + res.Tlr, 1.0; atol=1e-6)
    @test isapprox(res.Rll + res.Rrl + res.Tll + res.Trl, 1.0; atol=1e-6)
end

@testset "circular budget closes for a polarization-converting stack" begin
    # Half-wave plate at 45°: strong p<->s conversion in transmission. The old
    # scalar normalization N = (Tpp+Tss)/(|tpp|²+|tss|²) double-counted the
    # converted flux (Tpp/Tss were totals) and blew up when the co-polarized
    # amplitudes vanished at the half-wave point. The per-mode formulation
    # keeps every per-input circular budget exactly 1 for a lossless stack.
    no = λ -> 1.658 + 0.0im
    ne = λ -> 1.486 + 0.0im
    t_hw = 0.59 / (2 * (1.658 - 1.486))
    amb = Layer(λ -> 1.0 + 0.0im, 1.0)
    glass = Layer(λ -> 1.5 + 0.0im, 1.0)
    wp45 = Layer(no, no, ne, t_hw; euler=(π/4, π/2, 0.0))

    for (sub, θ) in ((amb, 0.0), (amb, 0.3), (glass, 0.3))
        res = transfer(0.59, [amb, wp45, sub]; θ=θ, basis=:circular)
        @test isapprox(res.Rrr + res.Rlr + res.Trr + res.Tlr, 1.0; atol=1e-10)
        @test isapprox(res.Rll + res.Rrl + res.Tll + res.Trl, 1.0; atol=1e-10)
        for v in (res.Trr, res.Tll, res.Trl, res.Tlr, res.Rrr, res.Rll, res.Rrl, res.Rlr)
            @test 0.0 <= v <= 1.0
        end
    end

    # A half-wave plate flips helicity in transmission: at the half-wave point
    # the co-handed Trr/Tll are small and the cross-handed Trl/Tlr carry the
    # transmitted power (linear-basis analogue: Tps ≈ 0.85).
    res = transfer(0.59, [amb, wp45, amb]; basis=:circular)
    @test res.Tlr > 0.8
    @test res.Trl > 0.8
    @test res.Trr < 0.05
    @test res.Tll < 0.05
end

@testset "circular TIR guard" begin
    # ξ = n_inc sinθ > n_substrate -> total internal reflection, no real flux.
    inc = Layer(λ -> 1.8 + 0im, 0.0)
    sub = Layer(λ -> 1.0 + 0im, 0.0)
    layers = [inc, sub]
    res = transfer(1.0, layers; θ=0.9, basis=:circular)

    for v in (res.Trr, res.Tll, res.Trl, res.Tlr)
        @test isapprox(v, 0.0; atol=1e-10)
        @test isfinite(v)
    end
    # All power reflected.
    @test isapprox(res.Rrr + res.Rlr, 1.0; atol=1e-6)
    @test isapprox(res.Rll + res.Rrl, 1.0; atol=1e-6)
end

@testset "sweep_angle basis=:circular" begin
    λs = [1.0, 1.2]
    θs = [0.0, 0.3]
    air = Layer(λs, fill(1.0, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(1.6, length(λs)), zeros(length(λs)), 0.25)
    layers = [air, film, air]

    spectra = sweep_angle(λs, θs, layers; basis=:circular)
    @test spectra isa CircularTransferResult{Matrix{Float64}}
    @test size(spectra.Rrl) == (length(θs), length(λs))
    @test size(spectra.Trr) == (length(θs), length(λs))

    # Each cell matches the scalar transfer call.
    for (i, θ) in enumerate(θs), (j, λ) in enumerate(λs)
        cell = transfer(λ, layers; θ=θ, basis=:circular)
        @test isapprox(spectra.Rrl[i, j], cell.Rrl; atol=1e-10)
        @test isapprox(spectra.Trr[i, j], cell.Trr; atol=1e-10)
    end

    # Default linear sweep unchanged.
    lin = sweep_angle(λs, θs, layers)
    @test lin isa TransferResult{Matrix{Float64}}
    @test lin.Rpp == sweep_angle(λs, θs, layers; basis=:linear).Rpp

    # Thread parity for the circular path.
    serial = sweep_angle(λs, θs, layers; basis=:circular, threads=false)
    @test spectra.Rrl == serial.Rrl
    @test spectra.Trr == serial.Trr
end

@testset "sweep_thickness basis=:circular" begin
    λs = [1.0, 1.1]
    air = Layer(λs, fill(1.0, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(1.5, length(λs)), zeros(length(λs)), 0.2)
    layers = [air, film, air]
    ts = [0.1, 0.2, 0.3]

    spectra = sweep_thickness(λs, ts, layers, 2; basis=:circular)
    @test spectra isa CircularTransferResult{Matrix{Float64}}
    @test size(spectra.Rrl) == (length(ts), length(λs))

    cell = transfer(λs[1], [air, Layer(λs, fill(1.5, length(λs)), zeros(length(λs)), ts[2]), air]; basis=:circular)
    @test isapprox(spectra.Rrl[2, 1], cell.Rrl; atol=1e-10)
end
