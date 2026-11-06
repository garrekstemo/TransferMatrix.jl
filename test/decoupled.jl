# 2×2 decoupled fast path: parity with the 4×4 :exp path on every output that
# consumes (M_sys, S), routing of :auto, and the :decoupled guard.

@testset "2×2 decoupled path" begin
    amb   = Layer(λ -> 1.0, 0.0)
    prism = Layer(λ -> 1.515, 0.0)
    glass = Layer(λ -> 1.52, 0.0)
    gold  = Layer(λ -> 0.18 + 3.4im, 0.05)
    tio2  = Layer(λ -> 2.4, 0.633 / 9.6)
    sio2  = Layer(λ -> 1.46, 0.633 / 5.84)
    lossy = Layer(λ -> 1.6 + 0.08im, 0.21)
    enz   = Layer(λ -> sqrt(complex(0.02 + 0.3im)), 0.12)
    sub_lossy = Layer(λ -> 3.5 + 0.2im, 0.0)
    uni0  = Layer(λ -> 1.55, λ -> 1.55, λ -> 1.72, 0.3)                 # optic axis ∥ z
    biax0 = Layer(λ -> 1.5, λ -> 1.62, λ -> 1.81, 0.27)                 # principal axes on lab axes
    mu_diag = Layer(λ -> 1.4, 0.2; mu=[1.3 0 0; 0 1.1 0; 0 0 0.9])
    mu_iso  = Layer(λ -> 1.4, 0.2; mu=1.7)
    dbr = vcat([amb], repeat([tio2, sio2], 8), [glass])
    graded = vcat([amb], [Layer(let n = 1.0 + 0.52 * (k - 0.5) / 60 + 0.001im; λ -> n; end, 0.3 / 60) for k in 1:60], [glass])

    stacks = [
        ("interface",      [amb, glass]),
        ("Kretschmann",    [prism, gold, amb]),
        ("DBR x8",         dbr),
        ("graded 60",      graded),
        ("lossy + ENZ",    [amb, lossy, enz, lossy, sub_lossy]),
        ("uniaxial, z axis", [amb, uni0, glass]),
        ("biaxial, lab axes", [amb, biax0, uni0, glass]),
        ("diagonal mu",    [amb, mu_diag, glass]),
        ("isotropic mu",   [amb, mu_iso, lossy, glass]),
    ]
    θs = (0.0, 0.3, 0.7, 1.1, 1.45)
    λs = (0.45, 0.633, 1.1)

    @testset "R/T parity with :exp — $name" for (name, ls) in stacks
        for λ in λs, θ in θs, μ in (1.0, 1.3)
            @test isdecoupled(λ, ls; θ=θ, μ=μ)
            a = transfer(λ, ls; θ=θ, μ=μ, method=:exp)
            b = transfer(λ, ls; θ=θ, μ=μ, method=:decoupled)
            for f in (:Tpp, :Tss, :Rpp, :Rss)
                @test isapprox(getfield(a, f), getfield(b, f); atol=1e-12)
            end
            for f in (:Tps, :Tsp, :Rps, :Rsp)
                @test getfield(b, f) == 0.0
                @test abs(getfield(a, f)) < 1e-12
            end
        end
    end

    @testset "amplitudes and system matrix parity — $name" for (name, ls) in stacks
        for λ in λs, θ in θs
            a = amplitude_coefficients(λ, ls; θ=θ, method=:exp)
            b = amplitude_coefficients(λ, ls; θ=θ, method=:decoupled)
            for f in (:rpp, :rss, :tpp, :tss)
                za, zb = getfield(a, f), getfield(b, f)
                @test abs(za - zb) <= 1e-10 * max(abs(za), 1.0)
            end
            Ma, _ = TransferMatrix._propagate(Val(:exp), λ, ls; θ=θ)
            Mb, _ = TransferMatrix._propagate(Val(:decoupled), λ, ls; θ=θ)
            @test maximum(abs, Ma - Mb) <= 1e-10 * max(maximum(abs, Ma), 1.0)
            @test all(iszero, (Mb[1,3], Mb[1,4], Mb[2,3], Mb[2,4], Mb[3,1], Mb[3,2], Mb[4,1], Mb[4,2]))
        end
    end

    @testset "complex θ beyond the ambient light line" begin
        for θ in (asin(1.2 + 0im), asin(1.05 + 0.01im))
            for f in (:rpp, :rss, :tpp, :tss)
                a = amplitude_coefficients(0.633, [prism, gold, amb]; θ=θ, method=:exp)
                b = amplitude_coefficients(0.633, [prism, gold, amb]; θ=θ, method=:decoupled)
                @test isapprox(getfield(a, f), getfield(b, f); rtol=1e-10)
            end
        end
    end

    @testset "circular basis parity" begin
        for θ in (0.0, 0.6)
            a = transfer(0.633, dbr; θ=θ, basis=:circular, method=:exp)
            b = transfer(0.633, dbr; θ=θ, basis=:circular, method=:decoupled)
            for f in fieldnames(typeof(a))
                @test isapprox(getfield(a, f), getfield(b, f); atol=1e-12)
            end
        end
    end

    @testset "diagonal sheets" begin
        graphene = Sheet(6.08e-5 + 1.0e-5im)
        tmdc = Sheet(; xx=1.0e-4 + 2.0e-5im, yy=8.0e-5 + 1.0e-5im)
        ls = [amb, sio2, lossy, glass]
        for sh in (Dict(1 => graphene), Dict(2 => tmdc), Dict(1 => tmdc, 3 => graphene)), θ in (0.0, 0.5, 1.2)
            @test isdecoupled(0.7, ls; θ=θ, sheets=sh)
            a = transfer(0.7, ls; θ=θ, sheets=sh, method=:exp)
            b = transfer(0.7, ls; θ=θ, sheets=sh, method=:decoupled)
            for f in (:Tpp, :Tss, :Rpp, :Rss)
                @test isapprox(getfield(a, f), getfield(b, f); atol=1e-12)
            end
        end
    end

    @testset "random isotropic stacks" begin
        # Deterministic low-discrepancy draws (no Random test dependency).
        u = let k = Ref(0); () -> (k[] += 1; mod(k[] * 0.6180339887498949, 1.0)); end
        for _ in 1:40
            N = 1 + floor(Int, 30u())
            ns = [complex(1.0 + 2.5u(), u() < 0.3 ? 0.5u() : 0.0) for _ in 1:N]
            ls = vcat([amb],
                      [Layer(let n = n; λ -> n; end, 0.4u()) for n in ns],
                      [Layer(let n = complex(1.3 + 2u(), 0.1u()); λ -> n; end, 0.0)])
            λ = 0.4 + u(); θ = 1.5u()
            a = transfer(λ, ls; θ=θ, method=:exp)
            b = transfer(λ, ls; θ=θ, method=:decoupled)
            for f in (:Tpp, :Tss, :Rpp, :Rss)
                @test isapprox(getfield(a, f), getfield(b, f); atol=1e-12)
            end
        end
    end

    @testset ":auto routing" begin
        rot   = Layer(λ -> 1.55, λ -> 1.55, λ -> 1.72, 0.3; euler=(π/4, π/4, 0.0))
        gyro  = Layer(λ -> 1.4, 0.2; mu=gyrotropic_tensor(2.0, 0.6))
        uniamb = Layer(λ -> 1.5, λ -> 1.5, λ -> 1.7, 0.0)
        hall = Sheet(; xx=1e-4, yy=1e-4, xy=2e-5, yx=-2e-5)
        coupled = [
            ("rotated uniaxial", [amb, rot, glass], nothing),
            ("gyrotropic mu",    [amb, gyro, glass], nothing),
            ("anisotropic ambient", [uniamb, sio2, glass], nothing),
            ("Hall sheet",       [amb, sio2, glass], Dict(1 => hall)),
        ]
        for (name, ls, sh) in coupled
            @test !isdecoupled(0.633, ls; θ=0.4, sheets=sh)
            @test_throws ArgumentError transfer(0.633, ls; θ=0.4, sheets=sh, method=:decoupled)
            a = transfer(0.633, ls; θ=0.4, sheets=sh, method=:exp)
            b = transfer(0.633, ls; θ=0.4, sheets=sh)            # :auto falls back to :exp
            for f in fieldnames(TransferResult)
                @test getfield(a, f) == getfield(b, f)
            end
        end
        # :auto on a decoupled stack is the 2×2 path itself
        a = transfer(0.633, dbr; θ=0.4, method=:decoupled)
        b = transfer(0.633, dbr; θ=0.4)
        for f in fieldnames(TransferResult)
            @test getfield(a, f) == getfield(b, f)
        end
    end

    @testset "sweeps default to :auto" begin
        λv = range(0.45, 0.9, length=7); θv = range(0.0, 1.2, length=5)
        a = sweep_angle(λv, θv, dbr; method=:exp, threads=false)
        b = sweep_angle(λv, θv, dbr; threads=false)
        @test maximum(abs, a.Rpp - b.Rpp) < 1e-12
        @test maximum(abs, a.Tss - b.Tss) < 1e-12
        ts = range(0.05, 0.2, length=4)
        a = sweep_thickness(λv, ts, dbr, 3; θ=0.3, method=:exp, threads=false)
        b = sweep_thickness(λv, ts, dbr, 3; θ=0.3, threads=false)
        @test maximum(abs, a.Rss - b.Rss) < 1e-12
        @test maximum(abs, a.Tpp - b.Tpp) < 1e-12
    end
end
