using BenchmarkTools
using TransferMatrix

BenchmarkTools.DEFAULT_PARAMETERS.samples = 5
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5.0

λ_0 = 1.0
λs = range(0.9, 1.1, length = 100)
θ = deg2rad(0.0)

n_air = 1.0
n_film = 1.5
n_sub = 1.45

d_film = λ_0 / (4 * n_film)
air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d_film)
sub = Layer(λs, fill(n_sub, length(λs)), zeros(length(λs)), 0.0)
layers = [air, film, sub]

println("Threads: ", Threads.nthreads())
println("calculate_tr single call")
trial = @benchmark calculate_tr($λ_0, $layers, $θ)
println("  min: ", BenchmarkTools.prettytime(minimum(trial).time),
        " | median: ", BenchmarkTools.prettytime(median(trial).time),
        " | allocs: ", trial.allocs,
        " | memory: ", Base.format_bytes(trial.memory))
