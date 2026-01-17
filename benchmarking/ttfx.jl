# Time-to-first-execution (TTFX) benchmark for PrecompileTools evaluation
#
# This script measures how long it takes from package load to first computation.
# Run from a fresh Julia session to get accurate TTFX measurements:
#
#   julia --project=benchmarking benchmarking/ttfx.jl
#
# For before/after comparison with PrecompileTools:
# 1. Run this script and record the times (before)
# 2. Add PrecompileTools with @compile_workload
# 3. Run `using TransferMatrix` once to trigger precompilation
# 4. Run this script again from a fresh session (after)

println("=" ^ 60)
println("TTFX Benchmark - TransferMatrix.jl")
println("=" ^ 60)
println()

# Measure package load time
t_load_start = time_ns()
using TransferMatrix
t_load_end = time_ns()
t_load = (t_load_end - t_load_start) / 1e9

println("Package load time: ", round(t_load, digits=3), " s")
println()

# Setup: create layers (uses simple dispersion functions, no RefractiveIndex.jl)
n_air = λ -> 1.0 + 0.0im
n_film = λ -> 1.5 + 0.0im
n_sub = λ -> 1.45 + 0.0im

λ_0 = 1.0  # μm
d_film = λ_0 / (4 * 1.5)  # quarter-wave thickness

air = Layer(n_air, 0.1)
film = Layer(n_film, d_film)
sub = Layer(n_sub, 0.5)
layers = [air, film, sub]

# Measure TTFX for transfer() - single point calculation
println("First call timings (TTFX):")
println("-" ^ 40)

t1_start = time_ns()
result1 = transfer(λ_0, layers)
t1_end = time_ns()
t_transfer = (t1_end - t1_start) / 1e9
println("transfer(λ, layers):         ", round(t_transfer, digits=4), " s")

# Measure TTFX for transfer() with angle
t2_start = time_ns()
result2 = transfer(λ_0, layers; θ=0.3)
t2_end = time_ns()
t_transfer_angle = (t2_end - t2_start) / 1e9
println("transfer(λ, layers; θ=0.3):  ", round(t_transfer_angle, digits=4), " s")

# Measure TTFX for efield()
t3_start = time_ns()
field = efield(λ_0, layers; dz=0.01)
t3_end = time_ns()
t_efield = (t3_end - t3_start) / 1e9
println("efield(λ, layers):           ", round(t_efield, digits=4), " s")

# Measure second call (should be fast - already compiled)
println()
println("Second call timings (compiled):")
println("-" ^ 40)

t4_start = time_ns()
result4 = transfer(λ_0, layers)
t4_end = time_ns()
t_transfer_2nd = (t4_end - t4_start) / 1e9
println("transfer(λ, layers):         ", round(t_transfer_2nd, digits=6), " s")

t5_start = time_ns()
field2 = efield(λ_0, layers; dz=0.01)
t5_end = time_ns()
t_efield_2nd = (t5_end - t5_start) / 1e9
println("efield(λ, layers):           ", round(t_efield_2nd, digits=6), " s")

println()
println("=" ^ 60)
println("Summary")
println("=" ^ 60)
println("Package load:              ", round(t_load, digits=3), " s")
println("First transfer() call:     ", round(t_transfer, digits=3), " s")
println("First efield() call:       ", round(t_efield, digits=3), " s")
println("Total TTFX (load + first): ", round(t_load + t_transfer, digits=3), " s")
println()
println("Speedup from compilation:")
println("  transfer(): ", round(t_transfer / t_transfer_2nd, digits=1), "x")
println("  efield():   ", round(t_efield / t_efield_2nd, digits=1), "x")
