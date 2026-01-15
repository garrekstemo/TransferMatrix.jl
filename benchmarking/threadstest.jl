using BenchmarkTools
using TransferMatrix

BenchmarkTools.DEFAULT_PARAMETERS.samples = 5
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5.0

# Simple air / film / substrate stack for a proof-of-concept benchmark.
λ_0 = 1.0
λs = range(0.9, 1.1, length = 100)
θs = range(0, 60, length = 100)

n_air = 1.0
n_film = 1.5
n_sub = 1.45

d_film = λ_0 / (4 * n_film)
d_air = 0.0
d_sub = 0.0

air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), d_air)
film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d_film)
sub = Layer(λs, fill(n_sub, length(λs)), zeros(length(λs)), d_sub)
layers = [air, film, sub]

println("Threads: ", Threads.nthreads())
println("sweep_angle grid: ", length(θs), " angles × ", length(λs), " wavelengths")
println("sweep_thickness grid: ", length(λs), " wavelengths × 100 thicknesses")

function run_bench(label, f)
    println(label)
    trial = @benchmark $f()
    min_time = minimum(trial).time
    med_time = median(trial).time
    println("  min: ", BenchmarkTools.prettytime(min_time),
            " | median: ", BenchmarkTools.prettytime(med_time),
            " | allocs: ", trial.allocs,
            " | memory: ", Base.format_bytes(trial.memory))
    return trial
end

trial_angle_threaded = run_bench(
    "sweep_angle threaded",
    () -> sweep_angle(λs, deg2rad.(θs), layers; threads=true),
)
trial_angle_serial = run_bench(
    "sweep_angle serial",
    () -> sweep_angle(λs, deg2rad.(θs), layers; threads=false),
)

ts = range(d_film * 0.8, d_film * 1.2, length = 100)
t_index = 2
trial_thick_threaded = run_bench(
    "sweep_thickness threaded",
    () -> sweep_thickness(λs, ts, layers, t_index, 0.0; threads=true),
)
trial_thick_serial = run_bench(
    "sweep_thickness serial",
    () -> sweep_thickness(λs, ts, layers, t_index, 0.0; threads=false),
)

angle_speedup = minimum(trial_angle_serial).time / minimum(trial_angle_threaded).time
thick_speedup = minimum(trial_thick_serial).time / minimum(trial_thick_threaded).time
println("Speedup (threaded vs serial): sweep_angle ×", round(angle_speedup; digits=2),
        ", sweep_thickness ×", round(thick_speedup; digits=2))
