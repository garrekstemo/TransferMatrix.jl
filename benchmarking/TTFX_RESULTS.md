# Time-to-First-Execution (TTFX) Benchmark Results

Benchmark comparing first-use latency before and after adding PrecompileTools.

## Test Configuration

- **Julia version**: 1.12
- **Date**: 2025-01-18
- **Hardware**: macOS (Darwin 25.2.0)
- **Benchmark script**: `benchmarking/ttfx.jl`

## Results

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| Package load | 0.91 s | 0.89 s | ~same |
| First `transfer()` call | 3.92 s | 0.022 s | **178x faster** |
| First `efield()` call | 0.53 s | 0.001 s | **530x faster** |
| **Total TTFX** | **4.83 s** | **0.92 s** | **5.3x faster** |

## Tradeoffs

- **Initial precompilation time**: Increased from ~2.3s to ~6s (one-time cost when installing/updating the package)
- **Package load time**: Unchanged (~0.9s)
- **Runtime performance**: Unchanged (already compiled code paths)

## Precompile Workload

The following operations are precompiled in `src/TransferMatrix.jl`:

```julia
@setup_workload begin
    n_air = λ -> 1.0 + 0.0im
    n_film = λ -> 1.5 + 0.0im
    n_sub = λ -> 1.45 + 0.0im

    λ_0 = 1.0
    d_film = λ_0 / 6.0

    @compile_workload begin
        air = Layer(n_air, 0.1)
        film = Layer(n_film, d_film)
        sub = Layer(n_sub, 0.5)
        layers = [air, film, sub]

        transfer(λ_0, layers)
        transfer(λ_0, layers; θ=0.3)
        efield(λ_0, layers; dz=0.01)
    end
end
```

## Running the Benchmark

To reproduce these results:

```bash
# Clear compiled cache (optional, for clean measurement)
rm -rf ~/.julia/compiled/v1.12/TransferMatrix*

# Trigger precompilation
julia --project -e 'using Pkg; Pkg.precompile()'

# Run TTFX benchmark (must be fresh Julia session)
julia --project=benchmarking benchmarking/ttfx.jl
```
