using BenchmarkTools
using TransferMatrix

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10
BenchmarkTools.DEFAULT_PARAMETERS.evals = 100

# Setup
λ = 1.0
n = 1.5 + 0.0im
layer = Layer(λ -> n, 0.1)
ξ = 0.0
μ = 1.0 + 0.0im

println("=" ^ 60)
println("Allocation breakdown for layer_matrices")
println("=" ^ 60)

# Benchmark layer_matrices (total)
println("\nlayer_matrices (total):")
trial = @benchmark TransferMatrix.layer_matrices($layer, $λ, $ξ, $μ)
println("  time: ", BenchmarkTools.prettytime(minimum(trial).time),
        " | allocs: ", trial.allocs,
        " | memory: ", Base.format_bytes(trial.memory))

# Setup intermediate values for individual function benchmarks
ε_i = TransferMatrix.dielectric_constant(n)
ε = TransferMatrix.dielectric_tensor(ε_i, ε_i, ε_i)
μ_tensor = TransferMatrix.permeability_tensor(μ, μ, μ)

println("\nIndividual functions:")
println("-" ^ 40)

# construct_M
print("construct_M:       ")
trial = @benchmark TransferMatrix.construct_M($ε, $μ_tensor)
println("allocs: ", lpad(trial.allocs, 3), " | ", BenchmarkTools.prettytime(minimum(trial).time))

M = TransferMatrix.construct_M(ε, μ_tensor)

# construct_a
print("construct_a:       ")
trial = @benchmark TransferMatrix.construct_a($ξ, $M)
println("allocs: ", lpad(trial.allocs, 3), " | ", BenchmarkTools.prettytime(minimum(trial).time))

a = TransferMatrix.construct_a(ξ, M)

# construct_Δ
print("construct_Δ:       ")
trial = @benchmark TransferMatrix.construct_Δ($ξ, $M, $a)
println("allocs: ", lpad(trial.allocs, 3), " | ", BenchmarkTools.prettytime(minimum(trial).time))

Δ = TransferMatrix.construct_Δ(ξ, M, a)

# calculate_q
print("calculate_q:       ")
trial = @benchmark TransferMatrix.calculate_q($Δ, $a)
println("allocs: ", lpad(trial.allocs, 3), " | ", BenchmarkTools.prettytime(minimum(trial).time))

q, S = TransferMatrix.calculate_q(Δ, a)

# calculate_γ
print("calculate_γ:       ")
trial = @benchmark TransferMatrix.calculate_γ($ξ, $q, $ε, $μ)
println("allocs: ", lpad(trial.allocs, 3), " | ", BenchmarkTools.prettytime(minimum(trial).time))

γ = TransferMatrix.calculate_γ(ξ, q, ε, μ)

# dynamical_matrix
print("dynamical_matrix:  ")
trial = @benchmark TransferMatrix.dynamical_matrix($ξ, $q, $γ, $μ)
println("allocs: ", lpad(trial.allocs, 3), " | ", BenchmarkTools.prettytime(minimum(trial).time))

# propagation_matrix
ω = 2π * 299792458 / λ
print("propagation_matrix:")
trial = @benchmark TransferMatrix.propagation_matrix($ω, $q)
println("allocs: ", lpad(trial.allocs, 3), " | ", BenchmarkTools.prettytime(minimum(trial).time))

println("\n" * "=" ^ 60)
