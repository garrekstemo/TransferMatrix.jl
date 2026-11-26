# ForwardDiff gradients through `transfer` and `amplitude_coefficients`, checked
# against central finite differences.

# Central-difference gradient of a scalar function f at p.
function _fd_gradient(f, p; h=1e-6)
    map(eachindex(p)) do i
        e = zeros(length(p)); e[i] = h * max(1.0, abs(p[i]))
        (f(p .+ e) - f(p .- e)) / (2e[i])
    end
end

@testset "ForwardDiff gradients" begin
    λ = 1e4 / 2100

    @testset "3-layer thicknesses, $m, θ=$θ, $fld" for m in (:auto, :decoupled, :exp),
            θ in (0.0, 0.35), fld in (:Tss, :Tpp, :Rss, :Rpp)
        layers(p) = [Layer(λ -> 1.0, 0.0), Layer(λ -> 2.4 + 0.01im, p[1]),
                     Layer(λ -> 1.35, p[2]), Layer(λ -> 1.5, p[3]), Layer(λ -> 1.0, 0.0)]
        f(p) = getproperty(transfer(λ, layers(p); θ=θ, method=m), fld)
        p = [0.5, 0.9, 0.45]
        @test f(p) isa Float64
        @test ForwardDiff.gradient(f, p) ≈ _fd_gradient(f, p) atol=1e-6
    end

    @testset "refractive indices, angle and wavelength" begin
        # n, k of a film and the substrate index
        layers(p) = [Layer(λ -> 1.0, 0.0), Layer(λ -> p[1] + p[2]*im, 0.7),
                     Layer(λ -> p[3], 0.0)]
        for m in (:auto, :exp), fld in (:Tss, :Tpp, :Rpp)
            f(p) = getproperty(transfer(λ, layers(p); θ=0.4, method=m), fld)
            p = [2.1, 0.05, 1.45]
            @test ForwardDiff.gradient(f, p) ≈ _fd_gradient(f, p) atol=1e-6
        end
        # θ and λ, with a dispersive (Lorentz) film
        stack = [Layer(λ -> 1.0, 0.0), Layer(lorentz(2.0, 1.2, 0.05; ε_inf=2.0), 0.3),
                 Layer(λ -> 1.5, 0.0)]
        for m in (:auto, :exp)
            gθ(θ) = transfer(1.0, stack; θ=θ, method=m).Tpp
            gλ(λ) = transfer(λ, stack; θ=0.3, method=m).Rss
            @test ForwardDiff.derivative(gθ, 0.3) ≈ _fd_gradient(x -> gθ(x[1]), [0.3])[1] atol=1e-6
            @test ForwardDiff.derivative(gλ, 1.0) ≈ _fd_gradient(x -> gλ(x[1]), [1.0])[1] atol=1e-6
        end
    end

    @testset "anisotropic, rotated and magnetic layers" begin
        no, ne = λ -> 1.66, λ -> 1.49
        # unrotated uniaxial film: 2×2 path; tilted in the plane of incidence: 2×2
        # via Δ; tilted out of it: 4×4 (s and p mix)
        for (ψtilt, m) in ((0.0, :auto), (0.6, :auto), (0.6, :exp))
            for φ in (0.0, 0.7)
                layers(p) = [Layer(λ -> 1.0, 0.0),
                             Layer(λ -> p[2], no, ne, p[1]; euler=(φ, p[3], 0.0)),
                             Layer(λ -> 1.5, 0.0)]
                for fld in (:Tpp, :Tss, :Tps, :Rsp)
                    f(p) = getproperty(transfer(λ, layers(p); θ=0.3, method=m), fld)
                    p = [0.8, 1.62, ψtilt]
                    @test ForwardDiff.gradient(f, p) ≈ _fd_gradient(f, p) atol=1e-6
                end
            end
        end
        # gyrotropic μ film (4×4 path)
        gyro(p) = [Layer(λ -> 1.0, 0.0),
                     Layer(λ -> 1.2, p[1]; mu=gyrotropic_tensor(p[2], p[3])),
                     Layer(λ -> 1.5, 0.0)]
        for fld in (:Tpp, :Tps, :Rss)
            f(p) = getproperty(transfer(λ, gyro(p); θ=0.3), fld)
            p = [0.4, 2.0, 0.6]
            @test ForwardDiff.gradient(f, p) ≈ _fd_gradient(f, p) atol=1e-6
        end
    end

    @testset "sheets and amplitude coefficients" begin
        stack(d) = [Layer(λ -> 1.0, 0.0), Layer(λ -> 1.8, d), Layer(λ -> 1.5, 0.0)]
        for m in (:auto, :exp)
            f(p) = transfer(λ, stack(p[1]); θ=0.2, method=m,
                            sheets=[1 => Sheet(p[2] * 1e-4 + 2e-4im)]).Tpp
            p = [0.6, 3.0]
            @test ForwardDiff.gradient(f, p) ≈ _fd_gradient(f, p) atol=1e-6
            g(p) = (c = amplitude_coefficients(λ, stack(p[1]); θ=0.2, method=m);
                    real(c.rpp) + imag(c.tss))
            @test ForwardDiff.gradient(g, [0.6]) ≈ _fd_gradient(g, [0.6]) atol=1e-6
            # dual index only: the first Float64 pass is redone with the dual type
            h(n) = real(amplitude_coefficients(λ, [Layer(λ -> 1.0, 0.0), Layer(λ -> n, 0.6),
                                                   Layer(λ -> 1.5, 0.0)]; θ=0.2, method=m).tpp)
            @test ForwardDiff.derivative(h, 1.8) ≈ _fd_gradient(x -> h(x[1]), [1.8])[1] atol=1e-6
        end
    end

    @testset "circular basis" begin
        layers(p) = [Layer(λ -> 1.0, 0.0), Layer(λ -> 1.66, λ -> 1.66, λ -> 1.49, p[1]; euler=(0.5, 0.9, 0.0)),
                     Layer(λ -> 1.5, 0.0)]
        f(p) = transfer(λ, layers(p); θ=0.2, basis=:circular).Trl
        @test ForwardDiff.gradient(f, [0.8]) ≈ _fd_gradient(f, [0.8]) atol=1e-6
    end

    @testset "Float64 results are unchanged" begin
        layers = [Layer(λ -> 1.0, 0.0), Layer(λ -> 2.4, 0.5), Layer(λ -> 1.5, 0.0)]
        @test transfer(λ, layers) isa TransferResult{Float64}
        @test transfer(λ, layers; method=:exp) isa TransferResult{Float64}
        @test transfer(λ, layers; method=:eig) isa TransferResult{Float64}
    end
end
