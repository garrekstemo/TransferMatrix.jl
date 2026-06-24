# Matrix-exponential propagator (issue #92): the public :exp backend, validated
# against the eigenmode path (:eig), the independent oracles in mu_reference.jl,
# energy conservation, and the interior-degeneracy robustness win.

if !isdefined(@__MODULE__, :MuReference)
    include("mu_reference.jl")
end
using .MuReference

@testset "matrix-exponential propagator (issue #92)" begin

    @testset "layer_transfer_exp == eig interior transfer (D-basis)" begin
        amb = Layer(λ -> 1.0, 1.0)
        cases = [
            ("isotropic",    Layer(λ -> 1.5, 0.3),                                      1.0),
            ("rot-uniaxial", Layer(λ -> 1.6, λ -> 1.6, λ -> 1.8, 0.4; euler=(0.4, 0.6, 0.2)), 1.0),
            ("scalar mu",    Layer(λ -> 1.5, 0.3),                                      2.0),
            ("tensor mu",    Layer(λ -> 1.5, 0.3; mu=gyrotropic_tensor(2.0, 0.6)),      1.0),
        ]
        for (_, L, μi) in cases, θ in (0.0, 0.6)
            λ = 1.0
            ω = 2π * TransferMatrix.c_0 / λ
            nx, _, _ = TransferMatrix.get_refractive_indices(amb, λ)
            ξ = sqrt(TransferMatrix.dielectric_constant(nx)) * sin(θ)
            D, P, _, _ = TransferMatrix.layer_matrices(L, λ, ξ, μi)
            Dm = SMatrix{4,4,ComplexF64}(D)
            T_eig = Dm * SMatrix{4,4,ComplexF64}(P(L.thickness)) * inv(Dm)
            T_exp = TransferMatrix.layer_transfer_exp(L, λ, ξ, ω, μi)
            @test maximum(abs, T_exp - T_eig) < 1e-12
        end
    end

end
