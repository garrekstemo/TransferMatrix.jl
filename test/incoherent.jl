# Incoherent (thick) layers: Katsidis & Siapkas intensity-matrix chaining.

# Independent reference: Macleod characteristic matrices (tilted admittances
# η_s = n cos θ, η_p = n / cos θ, T = Re(η_b)/Re(η_a)|t|²) for each coherent block,
# and the Katsidis & Siapkas 2×2 intensity-matrix product, written from the paper
# (Appl. Opt. 41, 3978, 2002) without package internals. `ns`, `ds` list the full
# stack (first and last semi-infinite), `inc` the incoherent layer indices,
# s = n₁ sin θ₁, pol ∈ (:s, :p). Returns (R, T).
function _ks_reference(λ, ns, ds, inc, s, pol)
    function cosθ(n)
        c = sqrt(complex(n^2 - s^2))
        (imag(c) < 0 || (imag(c) == 0 && real(c) < 0)) && (c = -c)
        return c / n
    end
    η(n) = pol === :s ? n * cosθ(n) : n / cosθ(n)
    function coherent(nb, db)   # (R, T) of the sub-stack nb, incident from nb[1]
        m = [1.0+0im 0; 0 1.0]
        for i in 2:length(nb)-1
            δ = 2π / λ * nb[i] * cosθ(nb[i]) * db[i]
            m = m * [cos(δ) -im*sin(δ)/η(nb[i]); -im*η(nb[i])*sin(δ) cos(δ)]
        end
        ηa, ηb = η(nb[1]), η(nb[end])
        B = m[1,1] + m[1,2] * ηb; Cc = m[2,1] + m[2,2] * ηb
        r = (ηa * B - Cc) / (ηa * B + Cc); t = 2ηa / (ηa * B + Cc)
        return abs2(r), real(ηb) / real(ηa) * abs2(t)
    end
    function intensity_matrix(ia, ib)
        Rf, Tf = coherent(ns[ia:ib], ds[ia:ib])
        Rb, Tb = coherent(reverse(ns[ia:ib]), reverse(ds[ia:ib]))
        return [1 -Rb; Rf (Tf * Tb - Rf * Rb)] ./ Tf
    end
    bounds = [1; inc; length(ns)]
    P = intensity_matrix(bounds[1], bounds[2])
    for j in 2:length(bounds)-1
        m = bounds[j]
        A = exp(-4π / λ * imag(ns[m] * cosθ(ns[m])) * ds[m])
        P = P * [1/A 0; 0 A] * intensity_matrix(m, bounds[j+1])
    end
    return P[2,1] / P[1,1], 1 / P[1,1]
end

@testset "incoherent layers" begin
    λ = 1e4 / 2000                     # 5 μm
    air = Layer(λ -> 1.0, 0.0)
    L(n, d) = Layer(λ -> n, d)

    @testset "matches the Katsidis–Siapkas matrix product" begin
        # film | 3 mm window | film | 1 mm window | film, both windows weakly absorbing
        ns = [1.0, 2.4, 1.35 + 0.002im, 1.40 + 2e-6im, 1.35, 2.4, 1.76 + 5e-6im, 2.4, 1.0]
        ds = [0.0, 0.5, 0.9, 3000.0, 0.7, 0.45, 1000.0, 0.6, 0.0]
        stack = [L(ns[i], ds[i]) for i in eachindex(ns)]
        for θ in (0.0, 0.35), (pol, Tf, Rf) in ((:s, :Tss, :Rss), (:p, :Tpp, :Rpp))
            res = transfer(λ, stack; θ=θ, incoherent=[4, 7])
            Rref, Tref = _ks_reference(λ, ns, ds, [4, 7], sin(θ), pol)
            @test getproperty(res, Tf) ≈ Tref rtol=1e-10
            @test getproperty(res, Rf) ≈ Rref rtol=1e-10
            @test res.Tps == res.Tsp == res.Rps == res.Rsp == 0
        end
    end

    @testset "equals the coherent spectrum averaged over a fringe" begin
        # lossless window: averaging the coherent T and R over one period of the
        # round-trip phase gives the incoherent result (to aliasing of order |r|^N)
        nC, θ = 1.40, 0.3
        stack(d) = [air, L(2.4, 0.5), L(1.35, 0.9), L(nC, d), L(1.35, 0.7), air]
        period = λ / (2 * sqrt(nC^2 - sin(θ)^2))
        inc = transfer(λ, stack(3000.0); θ=θ, incoherent=4)
        Nφ = 64
        coh = [transfer(λ, stack(3000.0 + k * period / Nφ); θ=θ) for k in 0:Nφ-1]
        for fld in (:Tss, :Tpp, :Rss, :Rpp)
            @test getproperty(inc, fld) ≈ sum(getproperty.(coh, fld)) / Nφ atol=1e-12
        end
        @test inc.Tss + inc.Rss ≈ 1 atol=1e-12
        @test inc.Tpp + inc.Rpp ≈ 1 atol=1e-12
    end

    @testset "limits" begin
        film = L(2.4 + 0.01im, 0.5)
        # a thick layer index-matched to the ambient adds nothing
        @test transfer(λ, [air, L(1.0, 5000.0), film, L(1.5, 0.0)]; θ=0.2, incoherent=2).Tpp ≈
              transfer(λ, [air, film, L(1.5, 0.0)]; θ=0.2).Tpp rtol=1e-12
        # an opaque thick layer: T → 0 and R → that of the front block, no NaN
        res = transfer(λ, [air, film, L(1.5 + 0.1im, 1e6), air]; incoherent=3)
        front = transfer(λ, [air, film, L(1.5 + 0.1im, 0.0)])
        @test res.Tss == 0
        @test res.Rss ≈ front.Rss rtol=1e-12
        # single incoherent index as an Int or in a tuple; empty means coherent
        st = [air, film, L(1.4, 3000.0), air]
        @test transfer(λ, st; incoherent=3).Tss == transfer(λ, st; incoherent=(3,)).Tss
        @test transfer(λ, st; incoherent=Int[]).Tss == transfer(λ, st).Tss
        # sheets inside a coherent block
        sh = [1 => Sheet(1e-3 + 2e-3im)]
        @test transfer(λ, [air, L(1.0, 3000.0), L(1.5, 0.0)]; sheets=[2 => Sheet(1e-3 + 2e-3im)], incoherent=2).Tpp ≈
              transfer(λ, [air, L(1.5, 0.0)]; sheets=sh).Tpp rtol=1e-12
        # sweeps pass it through
        sw = sweep_angle([λ, 1.1λ], [0.0, 0.2], st; incoherent=3, threads=false)
        @test sw.Tss[2, 1] == transfer(λ, st; θ=0.2, incoherent=3).Tss
    end

    @testset "argument checks" begin
        st = [air, L(2.4, 0.5), L(1.4, 3000.0), air]
        @test_throws ArgumentError transfer(λ, st; incoherent=1)
        @test_throws ArgumentError transfer(λ, st; incoherent=4)
        @test_throws ArgumentError transfer(λ, st; incoherent=3, method=:exp)
        @test_throws ArgumentError transfer(λ, st; incoherent=3, basis=:circular)
        aniso = [air, L(2.4, 0.5), Layer(λ -> 1.6, λ -> 1.6, λ -> 1.5, 3000.0), air]
        @test_throws ArgumentError transfer(λ, aniso; incoherent=3)
        coupled = [air, Layer(λ -> 1.6, λ -> 1.6, λ -> 1.5, 0.5; euler=(0.5, 0.7, 0.0)), L(1.4, 3000.0), air]
        @test_throws ArgumentError transfer(λ, coupled; θ=0.2, incoherent=3)
    end

    @testset "ForwardDiff through incoherent layers" begin
        stack(p) = [air, L(2.4, p[1]), L(1.35 + 0.002im, 0.9), Layer(λ -> p[3] + 1e-4im, p[2]),
                    L(1.35, 0.7), air]
        for fld in (:Tss, :Tpp, :Rpp)
            f(p) = getproperty(transfer(λ, stack(p); θ=0.2, incoherent=4), fld)
            p = [0.5, 3000.0, 1.40]
            @test ForwardDiff.gradient(f, p) ≈ _fd_gradient(f, p) atol=1e-6
        end
    end
end
