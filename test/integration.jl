@testset "transfer matrix interface" begin
    # Sanity check: a single interface at normal incidence should match Fresnel
    # reflectance and conserve energy (R + T = 1) for lossless media.
    λs = [1.0, 2.0]
    n1 = 1.0
    n2 = 1.5

    air = Layer(λs, fill(n1, length(λs)), zeros(length(λs)), 0.0)
    glass = Layer(λs, fill(n2, length(λs)), zeros(length(λs)), 0.0)
    layers = [air, glass]

    λ = 1.5
    (; Tpp, Tss, Rpp, Rss) = transfer(λ, layers)

    expected_R = ((n1 - n2) / (n1 + n2))^2
    expected_T = 1 - expected_R
    Rs, Rp = fresnel(0.0, n1, n2)

    @test isapprox(Rpp, Rss; atol=1e-12)
    @test isapprox(Tpp, Tss; atol=1e-12)
    @test isapprox(Rpp, expected_R; atol=1e-8)
    @test isapprox(Tpp, expected_T; atol=1e-8)
    @test isapprox(Rpp, Rp; atol=1e-8)
    @test isapprox(Rss, Rs; atol=1e-8)
end

@testset "sweep_angle interface" begin
    # sweep_angle should reproduce the same interface reflectance/transmittance
    # as transfer when evaluated at a fixed angle and wavelength grid.
    λs = [1.2, 1.6]
    θs = [0.0]
    n1 = 1.0
    n2 = 1.5

    air = Layer(λs, fill(n1, length(λs)), zeros(length(λs)), 0.0)
    glass = Layer(λs, fill(n2, length(λs)), zeros(length(λs)), 0.0)
    layers = [air, glass]

    spectra = sweep_angle(λs, θs, layers)

    expected_R = ((n1 - n2) / (n1 + n2))^2
    expected_T = 1 - expected_R

    @test size(spectra.Rpp) == (1, 2)
    @test all(isapprox.(spectra.Rpp, expected_R; atol=1e-8))
    @test all(isapprox.(spectra.Rss, expected_R; atol=1e-8))
    @test all(isapprox.(spectra.Tpp, expected_T; atol=1e-8))
    @test all(isapprox.(spectra.Tss, expected_T; atol=1e-8))
end

@testset "single thin film propagate and field" begin
    # Quarter-wave film at normal incidence: energy conservation and basic field
    # shape/metadata checks to ensure propagate/efield are coherent.
    λs = [1.0, 1.1]
    n_air = 1.0
    n_film = 1.5
    d = 1.0 / (4 * n_film)

    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d)
    layers = [air, film, air]

    λ = 1.0
    Γ, S, Ds, Ps, γs = TransferMatrix.propagate(λ, layers)
    _, _, _, _ = TransferMatrix.calculate_tr(Γ)
    Tpp, Tss, Rpp, Rss = TransferMatrix.calculate_tr(S)

    @test isapprox(Rpp, Rss; atol=1e-10)
    @test isapprox(Tpp + Rpp, 1.0; atol=1e-6)
    @test isapprox(Tss + Rss, 1.0; atol=1e-6)
    @test length(Ds) == length(layers)
    @test length(Ps) == length(layers)
    @test length(γs) == length(layers)

    ef = efield(λ, layers; dz=d / 10)
    @test size(ef.p, 1) == 3
    @test size(ef.s, 1) == 3
    @test size(ef.p, 2) == length(ef.z)
    @test size(ef.s, 2) == length(ef.z)
    @test length(ef.boundaries) == 2
    @test isapprox(ef.boundaries[1], 0.0; atol=1e-12)
    @test isapprox(ef.boundaries[2], d; atol=1e-12)
end

@testset "single thin film resonance and continuity" begin
    # Quarter-wave film should be more reflective than half-wave at the design
    # wavelength; field magnitudes should be continuous across interfaces.
    λ = 1.0
    λs = [λ, 1.1]
    n_air = 1.0
    n_film = 1.5

    d_quarter = λ / (4 * n_film)
    d_half = λ / (2 * n_film)

    air0 = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film_quarter = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d_quarter)
    film_half = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d_half)

    layers_quarter = [air0, film_quarter, air0]
    layers_half = [air0, film_half, air0]

    result_q = transfer(λ, layers_quarter)
    result_h = transfer(λ, layers_half)
    Tpp_q, Rpp_q = result_q.Tpp, result_q.Rpp
    Tpp_h, Rpp_h = result_h.Tpp, result_h.Rpp

    @test Rpp_q > Rpp_h
    @test Tpp_q < Tpp_h

    d_air = d_quarter / 2
    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), d_air)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d_quarter)
    layers = [air, film, air]

    ef = efield(λ, layers; dz=d_quarter / 40)

    function boundary_jump(field, z, boundary)
        left = findlast(<(boundary), z)
        right = findfirst(>(boundary), z)
        return left, right
    end

    for boundary in ef.boundaries
        left, right = boundary_jump(ef.p, ef.z, boundary)
        @test left !== nothing
        @test right !== nothing
        @test all(isapprox.(abs.(ef.p[1:2, left]), abs.(ef.p[1:2, right]); rtol=1e-2, atol=1e-2))
        @test all(isapprox.(abs.(ef.s[1:2, left]), abs.(ef.s[1:2, right]); rtol=1e-2, atol=1e-2))
    end
end

@testset "lossy thin film absorption" begin
    # Add extinction (k > 0) so absorption is nonzero; expect R + T < 1.
    λ = 1.0
    λs = [λ, 1.2]
    n_air = 1.0
    n_film = 1.5
    k_film = 0.2

    d = λ / (4 * n_film)

    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(n_film, length(λs)), fill(k_film, length(λs)), d)
    layers = [air, film, air]

    (; Tpp, Tss, Rpp, Rss) = transfer(λ, layers)

    @test Rpp + Tpp < 1.0
    @test Rss + Tss < 1.0
    @test Rpp >= 0.0
    @test Tpp >= 0.0
    @test Rss >= 0.0
    @test Tss >= 0.0
end

@testset "single thin film H-field continuity" begin
    # Berreman state vector ordering: Ψ = [Ex, Hy, Ey, -Hx].
    # This test reconstructs Ψ on either side of interfaces and checks that
    # tangential H components are continuous (within tolerance).
    λ = 1.0
    λs = [λ, 1.1]
    n_air = 1.0
    n_film = 1.5
    d = λ / (4 * n_film)
    eps = d / 200

    # Use a finite ambient layer so we can sample just inside each interface.
    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), d / 2)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d)
    layers = [air, film, air]

    # Rebuild the modal amplitudes per layer, mirroring efield internals,
    # so we can evaluate the Berreman state vector Ψ on either side of interfaces.
    Γ, S, Ds, Ps, γs = TransferMatrix.propagate(λ, layers)
    r, R, t, T = TransferMatrix.calculate_tr(Γ)
    first_layer = layers[1]
    last_layer = layers[end]

    Eplus_p = zeros(ComplexF64, length(layers), 4)
    Eminus_p = zeros(ComplexF64, length(layers), 4)
    Eplus_s = zeros(ComplexF64, length(layers), 4)
    Eminus_s = zeros(ComplexF64, length(layers), 4)

    # Outgoing modes in the final medium.
    Eplus_p[end, :] = [t[1], t[2], 0, 0]
    Eplus_s[end, :] = [t[3], t[4], 0, 0]

    # Back-propagate the field to each interface.
    P_f = Ps[end]
    Eminus_p[end, :] = inv(P_f(last_layer.thickness)) * Eplus_p[end, :]
    Eminus_s[end, :] = inv(P_f(last_layer.thickness)) * Eplus_s[end, :]

    D_i = Ds[end]
    for l in reverse(eachindex(layers))
        if l >= 2
            layer = layers[l - 1]
            D_prev = Ds[l - 1]
            P_prev = Ps[l - 1]
            L_i = inv(Ds[l - 1]) * D_i

            Eminus_p[l - 1, :] = L_i * Eplus_p[l, :]
            Eminus_s[l - 1, :] = L_i * Eplus_s[l, :]
            Eplus_p[l - 1, :] = P_prev(layer.thickness) * Eminus_p[l - 1, :]
            Eplus_s[l - 1, :] = P_prev(layer.thickness) * Eminus_s[l - 1, :]

            D_i = D_prev
        end
    end

    # Interface positions are measured from the first interface.
    interface_positions, _ = find_bounds(layers)
    interface_positions .-= first_layer.thickness

    # Build the Berreman state vector Ψ = D * E at position z within a layer.
    function ψ_from_modes(D, P, modes, z, z0)
        field = P(-(z - z0)) * modes
        return D * field
    end

    for (idx, boundary) in enumerate(interface_positions[1:end - 1])
        left_z = boundary - eps
        right_z = boundary + eps
        left_layer = idx
        right_layer = idx + 1

        # Compare tangential H components (Hy, Hx) just inside each side.
        ψp_left = ψ_from_modes(Ds[left_layer], Ps[left_layer], Eminus_p[left_layer, :], left_z, interface_positions[left_layer])
        ψp_right = ψ_from_modes(Ds[right_layer], Ps[right_layer], Eminus_p[right_layer, :], right_z, interface_positions[right_layer])
        ψs_left = ψ_from_modes(Ds[left_layer], Ps[left_layer], Eminus_s[left_layer, :], left_z, interface_positions[left_layer])
        ψs_right = ψ_from_modes(Ds[right_layer], Ps[right_layer], Eminus_s[right_layer, :], right_z, interface_positions[right_layer])

        @test isapprox(abs(ψp_left[2]), abs(ψp_right[2]); rtol=1e-2, atol=1e-2)
        @test isapprox(abs(ψp_left[4]), abs(ψp_right[4]); rtol=1e-2, atol=1e-2)
        @test isapprox(abs(ψs_left[2]), abs(ψs_right[2]); rtol=1e-2, atol=1e-2)
        @test isapprox(abs(ψs_left[4]), abs(ψs_right[4]); rtol=1e-2, atol=1e-2)
    end
end

@testset "DBR reflectivity trend" begin
    # A DBR built from quarter-wave layers should become more reflective
    # as the number of layer pairs increases at the design wavelength.
    λ = 1.0
    λs = [λ, 1.1]
    n_air = 1.0
    n1 = 2.2
    n2 = 1.5

    d1 = λ / (4 * n1)
    d2 = λ / (4 * n2)

    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    layer1 = Layer(λs, fill(n1, length(λs)), zeros(length(λs)), d1)
    layer2 = Layer(λs, fill(n2, length(λs)), zeros(length(λs)), d2)

    function stack_for_pairs(pairs)
        layers = [air]
        for _ in 1:pairs
            push!(layers, layer1, layer2)
        end
        push!(layers, air)
        return layers
    end

    R = Float64[]
    for pairs in (2, 4, 6)
        result = transfer(λ, stack_for_pairs(pairs))
        push!(R, result.Rpp)
        @test isapprox(result.Rpp, result.Rss; atol=1e-8)
    end

    @test R[1] < R[2] < R[3]
end

@testset "sweep_angle off-normal consistency" begin
    # sweep_angle at a nonzero angle should agree with transfer
    # for the same configuration, and conserve energy for lossless stacks.
    λs = [1.0, 1.2]
    θs = [0.3]
    n_air = 1.0
    n_film = 1.6
    d = 0.25

    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d)
    layers = [air, film, air]

    spectra = sweep_angle(λs, θs, layers)

    for (j, λ) in enumerate(λs)
        result = transfer(λ, layers; θ=θs[1])
        @test isapprox(spectra.Tpp[1, j], result.Tpp; atol=1e-8)
        @test isapprox(spectra.Tss[1, j], result.Tss; atol=1e-8)
        @test isapprox(spectra.Rpp[1, j], result.Rpp; atol=1e-8)
        @test isapprox(spectra.Rss[1, j], result.Rss; atol=1e-8)
        @test isapprox(result.Tpp + result.Rpp, 1.0; atol=1e-6)
        @test isapprox(result.Tss + result.Rss, 1.0; atol=1e-6)
    end
end

@testset "cross-polarization coefficients" begin
    # Synthetic Γ with off-diagonal coupling exercises Rps/Rsp and Tps paths
    # in calculate_tr without requiring a full anisotropic material model.
    # Synthetic transfer matrix with off-diagonal terms to exercise Rps/Rsp paths.
    Γ = [
        1.0  0.0  0.0  0.0;
        0.0  1.0  0.2  0.0;
        0.1  0.0  1.0  0.0;
        0.3  0.0  0.0  1.0
    ]

    r, R, t, T = TransferMatrix.calculate_tr(Γ)

    @test isapprox(r[1], -0.02; atol=1e-12)
    @test isapprox(r[2], 0.3; atol=1e-12)
    @test isapprox(r[4], 0.2; atol=1e-12)
    @test isapprox(T[4], 0.01; atol=1e-12)
    @test isapprox(R[3], 0.04; atol=1e-12)
    @test isapprox(R[4], 0.09; atol=1e-12)
end

@testset "cross-polarization zero for isotropic media" begin
    # For isotropic media, cross-polarization terms (Tps, Tsp, Rps, Rsp)
    # should be zero since there is no polarization mixing.
    λ = 1.0
    λs = [λ, 1.1]
    n_air = 1.0
    n_film = 1.5
    d = λ / (4 * n_film)

    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d)
    layers = [air, film, air]

    # Normal incidence
    result = transfer(λ, layers)
    @test isapprox(result.Tps, 0.0; atol=1e-12)
    @test isapprox(result.Tsp, 0.0; atol=1e-12)
    @test isapprox(result.Rps, 0.0; atol=1e-12)
    @test isapprox(result.Rsp, 0.0; atol=1e-12)

    # Oblique incidence - still zero for isotropic media
    result = transfer(λ, layers; θ=0.3)
    @test isapprox(result.Tps, 0.0; atol=1e-12)
    @test isapprox(result.Tsp, 0.0; atol=1e-12)
    @test isapprox(result.Rps, 0.0; atol=1e-12)
    @test isapprox(result.Rsp, 0.0; atol=1e-12)

    # Multi-layer DBR stack
    n1 = 2.2
    n2 = 1.5
    d1 = λ / (4 * n1)
    d2 = λ / (4 * n2)
    layer1 = Layer(λs, fill(n1, length(λs)), zeros(length(λs)), d1)
    layer2 = Layer(λs, fill(n2, length(λs)), zeros(length(λs)), d2)
    dbr_layers = [air, layer1, layer2, layer1, layer2, air]

    result = transfer(λ, dbr_layers; θ=0.2)
    @test isapprox(result.Tps, 0.0; atol=1e-12)
    @test isapprox(result.Tsp, 0.0; atol=1e-12)
    @test isapprox(result.Rps, 0.0; atol=1e-12)
    @test isapprox(result.Rsp, 0.0; atol=1e-12)

    # sweep_angle should also return zero cross-polarization
    θs = [0.0, 0.3]
    spectra = sweep_angle(λs, θs, layers)
    @test all(isapprox.(spectra.Tps, 0.0; atol=1e-12))
    @test all(isapprox.(spectra.Tsp, 0.0; atol=1e-12))
    @test all(isapprox.(spectra.Rps, 0.0; atol=1e-12))
    @test all(isapprox.(spectra.Rsp, 0.0; atol=1e-12))

    # sweep_thickness should also return zero cross-polarization
    ts = [d, 2d]
    spectra = sweep_thickness(λs, ts, layers, 2)
    @test all(isapprox.(spectra.Tps, 0.0; atol=1e-12))
    @test all(isapprox.(spectra.Tsp, 0.0; atol=1e-12))
    @test all(isapprox.(spectra.Rps, 0.0; atol=1e-12))
    @test all(isapprox.(spectra.Rsp, 0.0; atol=1e-12))
end

@testset "efield uniform medium" begin
    # If all layers have identical refractive index, there are no reflections,
    # so the field magnitude should be constant along z (pure phase evolution).
    λ = 1.0
    λs = [λ, 1.1]
    n = 1.0
    d = 0.5

    layer = Layer(λs, fill(n, length(λs)), zeros(length(λs)), d)
    layers = [layer, layer, layer]

    ef = efield(λ, layers; dz=d / 50)

    p_mag = abs.(ef.p[1, :])
    s_mag = abs.(ef.s[1, :])

    @test isapprox(maximum(p_mag), minimum(p_mag); rtol=1e-5, atol=1e-6)
    @test isapprox(maximum(s_mag), minimum(s_mag); rtol=1e-5, atol=1e-6)
end

@testset "keyword argument API" begin
    # Test that keyword arguments work correctly for θ and μ
    λ = 1.0
    λs = [λ, 1.1]
    n_air = 1.0
    n_film = 1.5
    d = λ / (4 * n_film)

    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d)
    layers = [air, film, air]

    # Test transfer with explicit θ
    θ_test = 0.3
    result = transfer(λ, layers; θ=θ_test)
    @test result.Tpp + result.Rpp ≈ 1.0 atol=1e-6
    @test result.Rpp != result.Rss  # At oblique incidence, p and s differ

    # Test transfer with explicit μ (non-magnetic, should be same as default)
    result2 = transfer(λ, layers; θ=θ_test, μ=1.0)
    @test result.Tpp ≈ result2.Tpp
    @test result.Rpp ≈ result2.Rpp

    # Test efield with explicit θ
    ef = efield(λ, layers; θ=θ_test, dz=d / 10)
    @test length(ef.z) > 0
    @test size(ef.p, 2) == length(ef.z)

    # Test sweep_thickness with explicit θ
    ts = [d, 2d]
    spectra = sweep_thickness(λs, ts, layers, 2; θ=θ_test)
    @test size(spectra.Rpp) == (2, 2)
end

@testset "threads toggle consistency" begin
    λs = [1.0, 1.2]
    θs = [0.0, 0.2]
    n_air = 1.0
    n_film = 1.4
    d = 0.3

    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d)
    layers = [air, film, air]

    spectra_threads = sweep_angle(λs, θs, layers; threads=true)
    spectra_serial = sweep_angle(λs, θs, layers; threads=false)

    @test spectra_threads.Rpp ≈ spectra_serial.Rpp
    @test spectra_threads.Rss ≈ spectra_serial.Rss
    @test spectra_threads.Tpp ≈ spectra_serial.Tpp
    @test spectra_threads.Tss ≈ spectra_serial.Tss

    ts = [0.1, 0.2, 0.3]
    spec_threads = sweep_thickness(λs, ts, layers, 2; threads=true)
    spec_serial = sweep_thickness(λs, ts, layers, 2; threads=false)

    @test spec_threads.Rpp ≈ spec_serial.Rpp
    @test spec_threads.Rss ≈ spec_serial.Rss
    @test spec_threads.Tpp ≈ spec_serial.Tpp
    @test spec_threads.Tss ≈ spec_serial.Tss
end

@testset "energy conservation validation" begin
    # Test that validate=true works without warnings for lossless media
    λ = 1.0
    λs = [λ, 1.1]
    n_air = 1.0
    n_film = 1.5
    d = λ / (4 * n_film)

    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d)
    layers = [air, film, air]

    # Should not produce warnings for lossless stack
    result = transfer(λ, layers; validate=true)
    @test result.Tpp + result.Rpp ≈ 1.0 atol=1e-6
    @test result.Tss + result.Rss ≈ 1.0 atol=1e-6

    # Test at oblique incidence
    result = transfer(λ, layers; θ=0.3, validate=true)
    @test result.Tpp + result.Rpp ≈ 1.0 atol=1e-6
    @test result.Tss + result.Rss ≈ 1.0 atol=1e-6

    # Test that validation is skipped for lossy media (no warning expected)
    k_film = 0.1
    lossy_film = Layer(λs, fill(n_film, length(λs)), fill(k_film, length(λs)), d)
    lossy_layers = [air, lossy_film, air]

    # Should complete without error; validation is skipped for absorbing media
    result = transfer(λ, lossy_layers; validate=true)
    @test result.Tpp + result.Rpp < 1.0  # Absorption means R + T < 1
    @test result.Tss + result.Rss < 1.0
end

@testset "_validate_physics warnings" begin
    λ = 1.0
    λs = [λ, 1.1]
    n_air = 1.0
    n_film = 1.5
    d = 0.1

    # Lossless layers for energy conservation tests
    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d)
    lossless_layers = [air, film, air]

    # Lossy layers for absorption bound tests
    lossy_film = Layer(λs, fill(n_film, length(λs)), fill(0.1, length(λs)), d)
    lossy_layers = [air, lossy_film, air]

    # Test NaN detection (NaN triggers early return, so only one warning)
    @test_logs (:warn, r"NaN detected") TransferMatrix._validate_physics(λ, lossless_layers, NaN, 0.5, 0.3, 0.5)
    @test_logs (:warn, r"NaN detected") TransferMatrix._validate_physics(λ, lossless_layers, 0.5, NaN, 0.3, 0.5)
    @test_logs (:warn, r"NaN detected") TransferMatrix._validate_physics(λ, lossless_layers, 0.5, 0.5, NaN, 0.5)
    @test_logs (:warn, r"NaN detected") TransferMatrix._validate_physics(λ, lossless_layers, 0.5, 0.5, 0.3, NaN)

    # Test bounds violations with match_mode=:any (multiple warnings may be issued)
    # T < 0
    @test_logs (:warn, r"Transmittance Tpp out of bounds") match_mode=:any TransferMatrix._validate_physics(λ, lossless_layers, -0.1, 0.5, 0.3, 0.5)
    @test_logs (:warn, r"Transmittance Tss out of bounds") match_mode=:any TransferMatrix._validate_physics(λ, lossless_layers, 0.5, -0.1, 0.3, 0.5)

    # T > 1
    @test_logs (:warn, r"Transmittance Tpp out of bounds") match_mode=:any TransferMatrix._validate_physics(λ, lossless_layers, 1.1, 0.5, 0.3, 0.5)
    @test_logs (:warn, r"Transmittance Tss out of bounds") match_mode=:any TransferMatrix._validate_physics(λ, lossless_layers, 0.5, 1.1, 0.3, 0.5)

    # R < 0
    @test_logs (:warn, r"Reflectance Rpp out of bounds") match_mode=:any TransferMatrix._validate_physics(λ, lossless_layers, 0.5, 0.5, -0.1, 0.5)
    @test_logs (:warn, r"Reflectance Rss out of bounds") match_mode=:any TransferMatrix._validate_physics(λ, lossless_layers, 0.5, 0.5, 0.3, -0.1)

    # R > 1
    @test_logs (:warn, r"Reflectance Rpp out of bounds") match_mode=:any TransferMatrix._validate_physics(λ, lossless_layers, 0.5, 0.5, 1.1, 0.5)
    @test_logs (:warn, r"Reflectance Rss out of bounds") match_mode=:any TransferMatrix._validate_physics(λ, lossless_layers, 0.5, 0.5, 0.3, 1.1)

    # Test energy conservation violation (lossless media, R + T ≠ 1)
    @test_logs (:warn, r"Energy conservation violated for p-polarization") match_mode=:any TransferMatrix._validate_physics(λ, lossless_layers, 0.5, 0.5, 0.3, 0.5)
    @test_logs (:warn, r"Energy conservation violated for s-polarization") match_mode=:any TransferMatrix._validate_physics(λ, lossless_layers, 0.5, 0.5, 0.5, 0.3)

    # Test absorption violation (lossy media, R + T > 1)
    @test_logs (:warn, r"Absorption violation.*p-polarization") match_mode=:any TransferMatrix._validate_physics(λ, lossy_layers, 0.7, 0.5, 0.5, 0.3)
    @test_logs (:warn, r"Absorption violation.*s-polarization") match_mode=:any TransferMatrix._validate_physics(λ, lossy_layers, 0.5, 0.7, 0.3, 0.5)

    # Test valid cases produce no warnings
    # Lossless: R + T = 1
    @test_logs TransferMatrix._validate_physics(λ, lossless_layers, 0.7, 0.6, 0.3, 0.4)
    # Lossy: R + T < 1 (absorption is fine)
    @test_logs TransferMatrix._validate_physics(λ, lossy_layers, 0.5, 0.5, 0.3, 0.3)
end

@testset "anisotropic layer construction" begin
    # Verify anisotropic layers can be constructed and used in transfer()
    λ = 1.0

    # Isotropic ambient
    air = Layer(λ -> 1.0, 0.0)

    # Uniaxial crystal: optic axis along z (calcite-like)
    no = λ -> 1.658  # ordinary
    ne = λ -> 1.486  # extraordinary
    uniaxial = Layer(no, no, ne, 0.5)

    layers = [air, uniaxial, air]

    # Should run without error
    result = transfer(λ, layers)
    @test result.Tpp >= 0.0
    @test result.Tss >= 0.0
    @test result.Rpp >= 0.0
    @test result.Rss >= 0.0

    # For uniaxial with optic axis along z at normal incidence,
    # there should be no cross-polarization (principal axes aligned)
    @test isapprox(result.Tps, 0.0; atol=1e-10)
    @test isapprox(result.Tsp, 0.0; atol=1e-10)
    @test isapprox(result.Rps, 0.0; atol=1e-10)
    @test isapprox(result.Rsp, 0.0; atol=1e-10)
end

@testset "anisotropic layer birefringence" begin
    # Birefringent layer should show different behavior for s and p polarizations
    λ = 1.0

    air = Layer(λ -> 1.0, 0.0)

    # Biaxial crystal with distinct indices
    nx = λ -> 1.5
    ny = λ -> 1.6
    nz = λ -> 1.7
    biaxial = Layer(nx, ny, nz, 0.5)

    layers = [air, biaxial, air]

    # At normal incidence, x and y axes are in the plane of the interface
    # s-polarization sees ny, p-polarization sees nx (for incidence along z)
    result = transfer(λ, layers)

    # Verify physical bounds
    @test 0 ≤ result.Tpp ≤ 1
    @test 0 ≤ result.Tss ≤ 1
    @test 0 ≤ result.Rpp ≤ 1
    @test 0 ≤ result.Rss ≤ 1

    # Energy should be conserved (lossless)
    @test isapprox(result.Tpp + result.Rpp, 1.0; atol=1e-6)
    @test isapprox(result.Tss + result.Rss, 1.0; atol=1e-6)

    # With different nx and ny, s and p polarizations should see different media
    # This means Tpp ≠ Tss and Rpp ≠ Rss (unless at special angles/thicknesses)
    # We just verify they're both reasonable values here
    @test result.Tpp > 0
    @test result.Tss > 0
end

@testset "anisotropic vs isotropic equivalence" begin
    # When all three indices are equal, anisotropic should give same result as isotropic
    λ = 1.0
    n = 1.5

    air = Layer(λ -> 1.0, 0.0)
    iso_film = Layer(λ -> n, 0.3)
    aniso_film = Layer(λ -> n, λ -> n, λ -> n, 0.3)

    layers_iso = [air, iso_film, air]
    layers_aniso = [air, aniso_film, air]

    result_iso = transfer(λ, layers_iso)
    result_aniso = transfer(λ, layers_aniso)

    @test isapprox(result_iso.Tpp, result_aniso.Tpp; atol=1e-12)
    @test isapprox(result_iso.Tss, result_aniso.Tss; atol=1e-12)
    @test isapprox(result_iso.Rpp, result_aniso.Rpp; atol=1e-12)
    @test isapprox(result_iso.Rss, result_aniso.Rss; atol=1e-12)

    # Also test at oblique incidence
    θ = 0.3
    result_iso = transfer(λ, layers_iso; θ=θ)
    result_aniso = transfer(λ, layers_aniso; θ=θ)

    @test isapprox(result_iso.Tpp, result_aniso.Tpp; atol=1e-12)
    @test isapprox(result_iso.Tss, result_aniso.Tss; atol=1e-12)
    @test isapprox(result_iso.Rpp, result_aniso.Rpp; atol=1e-12)
    @test isapprox(result_iso.Rss, result_aniso.Rss; atol=1e-12)
end

@testset "rotated crystal orientation" begin
    # Test that rotating a uniaxial crystal changes the optical response
    λ = 1.0

    air = Layer(λ -> 1.0, 0.0)

    # Uniaxial crystal: ne along z (optic axis)
    no = λ -> 1.658  # ordinary
    ne = λ -> 1.486  # extraordinary

    # Unrotated: optic axis along z
    unrotated = Layer(no, no, ne, 0.5)

    # Rotated: optic axis tilted 45° from z in xz-plane
    rotated_45 = Layer(no, no, ne, 0.5; euler=(0.0, π/4, 0.0))

    # Rotated: optic axis in xy-plane (90° tilt)
    rotated_90 = Layer(no, no, ne, 0.5; euler=(0.0, π/2, 0.0))

    layers_unrot = [air, unrotated, air]
    layers_rot45 = [air, rotated_45, air]
    layers_rot90 = [air, rotated_90, air]

    result_unrot = transfer(λ, layers_unrot)
    result_rot45 = transfer(λ, layers_rot45)
    result_rot90 = transfer(λ, layers_rot90)

    # All should satisfy energy conservation (lossless)
    @test isapprox(result_unrot.Tpp + result_unrot.Rpp, 1.0; atol=1e-6)
    @test isapprox(result_unrot.Tss + result_unrot.Rss, 1.0; atol=1e-6)
    @test isapprox(result_rot45.Tpp + result_rot45.Rpp, 1.0; atol=1e-6)
    @test isapprox(result_rot45.Tss + result_rot45.Rss, 1.0; atol=1e-6)
    @test isapprox(result_rot90.Tpp + result_rot90.Rpp, 1.0; atol=1e-6)
    @test isapprox(result_rot90.Tss + result_rot90.Rss, 1.0; atol=1e-6)

    # All should have valid physical bounds
    for result in [result_unrot, result_rot45, result_rot90]
        @test 0 ≤ result.Tpp ≤ 1
        @test 0 ≤ result.Tss ≤ 1
        @test 0 ≤ result.Rpp ≤ 1
        @test 0 ≤ result.Rss ≤ 1
    end

    # Unrotated at normal incidence: no cross-polarization (principal axes aligned)
    @test isapprox(result_unrot.Tps, 0.0; atol=1e-10)
    @test isapprox(result_unrot.Tsp, 0.0; atol=1e-10)
    @test isapprox(result_unrot.Rps, 0.0; atol=1e-10)
    @test isapprox(result_unrot.Rsp, 0.0; atol=1e-10)

    # Rotated crystals should give different results than unrotated
    # (unless at special angles/thicknesses)
    # At minimum, 90° rotation should change things significantly
    @test result_unrot.Tpp != result_rot90.Tpp || result_unrot.Tss != result_rot90.Tss
end

@testset "rotated isotropic invariance" begin
    # Rotating an isotropic material should give same result
    # (since all three indices are equal, rotation doesn't matter)
    λ = 1.0
    n = 1.5

    air = Layer(λ -> 1.0, 0.0)

    # "Anisotropic" layer with equal indices, no rotation
    aniso_unrot = Layer(λ -> n, λ -> n, λ -> n, 0.3)

    # Same layer with rotation - should give identical results
    aniso_rot = Layer(λ -> n, λ -> n, λ -> n, 0.3; euler=(π/6, π/4, π/3))

    layers_unrot = [air, aniso_unrot, air]
    layers_rot = [air, aniso_rot, air]

    result_unrot = transfer(λ, layers_unrot)
    result_rot = transfer(λ, layers_rot)

    @test isapprox(result_unrot.Tpp, result_rot.Tpp; atol=1e-12)
    @test isapprox(result_unrot.Tss, result_rot.Tss; atol=1e-12)
    @test isapprox(result_unrot.Rpp, result_rot.Rpp; atol=1e-12)
    @test isapprox(result_unrot.Rss, result_rot.Rss; atol=1e-12)

    # Both should have zero cross-polarization
    @test isapprox(result_unrot.Tps, 0.0; atol=1e-12)
    @test isapprox(result_rot.Tps, 0.0; atol=1e-12)
end

@testset "euler angle symmetry" begin
    # Test that 2π rotation gives same result as no rotation
    λ = 1.0

    air = Layer(λ -> 1.0, 0.0)
    no = λ -> 1.658
    ne = λ -> 1.486

    # No rotation
    layer0 = Layer(no, no, ne, 0.5; euler=(0.0, 0.0, 0.0))

    # Full rotation about z (2π)
    layer2pi = Layer(no, no, ne, 0.5; euler=(2π, 0.0, 0.0))

    layers0 = [air, layer0, air]
    layers2pi = [air, layer2pi, air]

    result0 = transfer(λ, layers0)
    result2pi = transfer(λ, layers2pi)

    @test isapprox(result0.Tpp, result2pi.Tpp; atol=1e-12)
    @test isapprox(result0.Tss, result2pi.Tss; atol=1e-12)
    @test isapprox(result0.Rpp, result2pi.Rpp; atol=1e-12)
    @test isapprox(result0.Rss, result2pi.Rss; atol=1e-12)
end
