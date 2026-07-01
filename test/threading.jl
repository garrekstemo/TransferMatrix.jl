# Regression tests for thread safety of the sweep APIs.
#
# sweep_angle / sweep_thickness use Threads.@threads over a shared `layers`
# vector, so every layer's dispersion function is evaluated concurrently.
# Tabulated-data layers built by `refractive_index(λs, ns, ks)` previously
# captured DataInterpolations interpolants whose Guesser wrote a search-hint
# RefValue on every evaluation — a formal data race across threads. These tests
# put a tabulated layer in threaded sweeps and pin the results against serial
# per-point `transfer` calls. Run the suite with 2+ threads (julia -t 2) for
# the concurrency to be exercised.

@testset "threaded sweeps with tabulated-data layer" begin
    if Threads.nthreads() == 1
        @info "threading.jl: running with 1 thread; start Julia with -t 2+ to exercise concurrency"
    end

    # Non-monotonic n(λ) so successive interpolation lookups jump between knots
    λ_data = collect(range(0.4, 1.6, length=25))
    n_data = 1.4 .+ 0.1 .* sin.(2π .* λ_data)
    k_data = 0.05 .+ 0.02 .* cos.(2π .* λ_data)

    air = Layer(λ -> 1.0, 0.1)
    film = Layer(λ_data, n_data, k_data, 0.3)
    sub = Layer(λ -> 1.45, 0.1)
    layers = [air, film, sub]

    λs = collect(range(0.5, 1.5, length=40))
    θs = collect(range(0.0, 1.2, length=16))

    result = sweep_angle(λs, θs, layers)

    # Serial reference computed one (θ, λ) point at a time on this thread
    ref_Tss = [transfer(λ, layers; θ=θ).Tss for θ in θs, λ in λs]
    ref_Tpp = [transfer(λ, layers; θ=θ).Tpp for θ in θs, λ in λs]
    ref_Rss = [transfer(λ, layers; θ=θ).Rss for θ in θs, λ in λs]
    ref_Rpp = [transfer(λ, layers; θ=θ).Rpp for θ in θs, λ in λs]

    @test result.Tss == ref_Tss
    @test result.Tpp == ref_Tpp
    @test result.Rss == ref_Rss
    @test result.Rpp == ref_Rpp

    # sweep_thickness threads over the same shared tabulated layer
    ts = collect(range(0.1, 0.5, length=12))
    result_t = sweep_thickness(λs, ts, layers, 2; θ=0.2)
    ref_t = [transfer(λ, [air, Layer(λ_data, n_data, k_data, t), sub]; θ=0.2).Tss
             for t in ts, λ in λs]
    @test result_t.Tss == ref_t
end
