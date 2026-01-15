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
    Tpp, Tss, Rpp, Rss = calculate_tr(λ, layers, 0.0)

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

@testset "angle_resolved interface" begin
    # angle_resolved should reproduce the same interface reflectance/transmittance
    # as calculate_tr when evaluated at a fixed angle and wavelength grid.
    λs = [1.2, 1.6]
    θs = [0.0]
    n1 = 1.0
    n2 = 1.5

    air = Layer(λs, fill(n1, length(λs)), zeros(length(λs)), 0.0)
    glass = Layer(λs, fill(n2, length(λs)), zeros(length(λs)), 0.0)
    layers = [air, glass]

    spectra = angle_resolved(λs, θs, layers)

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
    # shape/metadata checks to ensure propagate/electric_field are coherent.
    λs = [1.0, 1.1]
    n_air = 1.0
    n_film = 1.5
    d = 1.0 / (4 * n_film)

    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d)
    layers = [air, film, air]

    λ = 1.0
    Γ, S, Ds, Ps, γs = TransferMatrix.propagate(λ, layers, 0.0, 1.0 + 0.0im)
    _, _, _, _ = calculate_tr(Γ)
    Tpp, Tss, Rpp, Rss = calculate_tr(S)

    @test isapprox(Rpp, Rss; atol=1e-10)
    @test isapprox(Tpp + Rpp, 1.0; atol=1e-6)
    @test isapprox(Tss + Rss, 1.0; atol=1e-6)
    @test length(Ds) == length(layers)
    @test length(Ps) == length(layers)
    @test length(γs) == length(layers)

    ef = electric_field(λ, layers, 0.0; dz=d / 10)
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

    Tpp_q, _, Rpp_q, _ = calculate_tr(λ, layers_quarter, 0.0)
    Tpp_h, _, Rpp_h, _ = calculate_tr(λ, layers_half, 0.0)

    @test Rpp_q > Rpp_h
    @test Tpp_q < Tpp_h

    d_air = d_quarter / 2
    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), d_air)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d_quarter)
    layers = [air, film, air]

    ef = electric_field(λ, layers, 0.0; dz=d_quarter / 40)

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

    Tpp, Tss, Rpp, Rss = calculate_tr(λ, layers, 0.0)

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

    # Rebuild the modal amplitudes per layer, mirroring electric_field internals,
    # so we can evaluate the Berreman state vector Ψ on either side of interfaces.
    Γ, S, Ds, Ps, γs = TransferMatrix.propagate(λ, layers, 0.0, 1.0 + 0.0im)
    r, R, t, T = calculate_tr(Γ)
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
        Tpp, Tss, Rpp, Rss = calculate_tr(λ, stack_for_pairs(pairs), 0.0)
        push!(R, Rpp)
        @test isapprox(Rpp, Rss; atol=1e-8)
    end

    @test R[1] < R[2] < R[3]
end

@testset "angle_resolved off-normal consistency" begin
    # angle_resolved at a nonzero angle should agree with calculate_tr
    # for the same configuration, and conserve energy for lossless stacks.
    λs = [1.0, 1.2]
    θs = [0.3]
    n_air = 1.0
    n_film = 1.6
    d = 0.25

    air = Layer(λs, fill(n_air, length(λs)), zeros(length(λs)), 0.0)
    film = Layer(λs, fill(n_film, length(λs)), zeros(length(λs)), d)
    layers = [air, film, air]

    spectra = angle_resolved(λs, θs, layers)

    for (j, λ) in enumerate(λs)
        Tpp, Tss, Rpp, Rss = calculate_tr(λ, layers, θs[1])
        @test isapprox(spectra.Tpp[1, j], Tpp; atol=1e-8)
        @test isapprox(spectra.Tss[1, j], Tss; atol=1e-8)
        @test isapprox(spectra.Rpp[1, j], Rpp; atol=1e-8)
        @test isapprox(spectra.Rss[1, j], Rss; atol=1e-8)
        @test isapprox(Tpp + Rpp, 1.0; atol=1e-6)
        @test isapprox(Tss + Rss, 1.0; atol=1e-6)
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

    r, R, t, T = calculate_tr(Γ)

    @test isapprox(r[1], -0.02; atol=1e-12)
    @test isapprox(r[2], 0.3; atol=1e-12)
    @test isapprox(r[4], 0.2; atol=1e-12)
    @test isapprox(T[4], 0.01; atol=1e-12)
    @test isapprox(R[3], 0.04; atol=1e-12)
    @test isapprox(R[4], 0.09; atol=1e-12)
end

@testset "electric_field uniform medium" begin
    # If all layers have identical refractive index, there are no reflections,
    # so the field magnitude should be constant along z (pure phase evolution).
    λ = 1.0
    λs = [λ, 1.1]
    n = 1.0
    d = 0.5

    layer = Layer(λs, fill(n, length(λs)), zeros(length(λs)), d)
    layers = [layer, layer, layer]

    ef = electric_field(λ, layers, 0.0; dz=d / 50)

    p_mag = abs.(ef.p[1, :])
    s_mag = abs.(ef.s[1, :])

    @test isapprox(maximum(p_mag), minimum(p_mag); rtol=1e-5, atol=1e-6)
    @test isapprox(maximum(s_mag), minimum(s_mag); rtol=1e-5, atol=1e-6)
end
