@testset "find_bounds" begin
    l1 = TransferMatrix.Layer(RefractiveMaterial("main", "Au", "Rakic-LD"), 1)
    l2 = TransferMatrix.Layer(RefractiveMaterial("main", "SiO2", "Malitson"), 2)
    l3 = TransferMatrix.Layer(RefractiveMaterial("main", "Au", "Rakic-LD"), 3)

    layers = [l1, l2, l3]

    interface_positions, total_thickness = TransferMatrix.find_bounds(layers)

    @test interface_positions == [1.0, 3.0, 6.0]
    @test total_thickness == 6.0

    single_layer = [TransferMatrix.Layer(RefractiveMaterial("main", "SiO2", "Malitson"), 1.5)]
    single_positions, single_total = TransferMatrix.find_bounds(single_layer)

    @test single_positions == [1.5]
    @test single_total == 1.5
end


# Recursively check whether a Base.RefValue is reachable from `x`. Pins the
# refractive_index thread-safety fix: a RefValue captured by a dispersion closure
# means per-evaluation mutable state, which the threaded sweeps race on.
function captures_refvalue(x, seen=Base.IdSet{Any}())
    x isa Base.RefValue && return true
    T = typeof(x)
    (isbitstype(T) || x isa Union{Symbol, String, Module, Type}) && return false
    x in seen && return false
    push!(seen, x)
    if x isa AbstractArray
        isbitstype(eltype(x)) && return false
        return any(i -> isassigned(x, i) && captures_refvalue(x[i], seen), eachindex(x))
    end
    isstructtype(T) || return false
    return any(i -> isdefined(x, i) && captures_refvalue(getfield(x, i), seen), 1:fieldcount(T))
end

@testset "refractive_index" begin
    air = RefractiveMaterial("other", "air", "Ciddor")
    au = RefractiveMaterial("main", "Au", "Rakic-LD")

    @test_throws ArgumentError extinction(air, 1.0)
    @test TransferMatrix.refractive_index(air)(1.0) == dispersion(air, 1.0) + 0.0im
    @test TransferMatrix.refractive_index(air)(1.0) == 1.0002741661312147 + 0.0im
    @test TransferMatrix.refractive_index(au)(1.0) == 0.2557301597051597 + 5.986408108108109im

    λs = [1.0, 2.0, 3.0]
    ns = [1.5, 2.0, 2.5]
    ks = [0.5, 1.0, 1.5]
    refractive_index_func = refractive_index(λs, ns, ks)

    @test refractive_index_func(1.0) == 1.5 + 0.5im
    @test refractive_index_func(2.0) == 2.0 + 1.0im
    @test refractive_index_func(3.0) == 2.5 + 1.5im

    # Linear interpolation between knots
    @test refractive_index_func(1.5) ≈ 1.75 + 0.75im
    @test refractive_index_func(2.25) ≈ 2.125 + 1.125im

    # Out-of-domain wavelengths throw rather than silently extrapolating
    @test_throws DomainError refractive_index_func(0.5)
    @test_throws DomainError refractive_index_func(3.5)

    # Construction-time validation
    @test_throws ArgumentError refractive_index([3.0, 2.0, 1.0], ns, ks)  # unsorted λs
    @test_throws ArgumentError refractive_index(λs, [1.5, 2.0], ks)      # length mismatch
    @test_throws ArgumentError refractive_index([1.0], [1.5], [0.5])     # too few points
    # Duplicate wavelengths are data errors: a duplicated final knot made the
    # endpoint evaluate to NaN (zero-width segment), so reject at construction.
    @test_throws ArgumentError refractive_index([1.0, 2.0, 3.0, 3.0], [1.5, 2.0, 2.5, 9.0], [0.0, 0.0, 0.0, 0.0])

    # Thread-safety regression: the closure must not capture mutable evaluation
    # state. DataInterpolations' LinearInterpolation wrote a Guesser idx_prev
    # RefValue on every call — a data race in the threaded sweeps.
    @test !captures_refvalue(refractive_index_func)
end


@testset "dielectric_constant" begin
    au = RefractiveMaterial("main", "Au", "Rakic-LD")
    l = TransferMatrix.Layer(au, 1.0)

    ε = TransferMatrix.dielectric_constant.([1.0, 1.0, 1.0], [2.0, 2.0, 2.0])
    @test real(ε) == [-3.0, -3.0, -3.0]
    @test imag(ε) == [4.0, 4.0, 4.0]

    @test TransferMatrix.dielectric_constant(1.0, 2.0) == -3.0 + 4.0im
    @test TransferMatrix.dielectric_constant(1.0 + 2.0im) == -3.0 + 4.0im

end


@testset "dielectric_tensor" begin
    @test TransferMatrix.dielectric_tensor(1.0, 1.0, 1.0) == [1.0 0 0; 0 1.0 0; 0 0 1.0]
    @test TransferMatrix.dielectric_tensor(1.0 + 1.0im, 1.0 + 1.0im, 1.0) == [complex(1.0, 1.0) 0 0; 0 complex(1.0, 1.0) 0; 0 0 1.0]
end


@testset "permeability_tensor" begin
    μ1, μ2, μ3 = 1.0 + 1.0im, 2.0 + 2.0im, 3.0 + 3.0im
    expected_tensor = Diagonal(SVector{3, ComplexF64}(μ1, μ2, μ3))

    @test TransferMatrix.permeability_tensor(μ1, μ2, μ3) == expected_tensor
end


@testset "construct_constitutive" begin
    ε_i = (1.0 + 2.0im)^2
    μ_i = 1.0 + 0.0im
    ε = Diagonal([ε_i, ε_i, ε_i])
    μ = Diagonal(fill(μ_i, 3))
    ρ1 = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    ρ2 = [9.0 8.0 7.0; 6.0 5.0 4.0; 3.0 2.0 1.0]

    M1 = TransferMatrix.construct_constitutive(ε, μ, ρ1, ρ2)
    M1_true = [
        ε_i 0.0 0.0 1.0 2.0 3.0;
        0.0 ε_i 0.0 4.0 5.0 6.0;
        0.0 0.0 ε_i 7.0 8.0 9.0;
        9.0 8.0 7.0 μ_i 0.0 0.0;
        6.0 5.0 4.0 0.0 μ_i 0.0;
        3.0 2.0 1.0 0.0 0.0 μ_i
    ]

    M2 = TransferMatrix.construct_constitutive(ε, μ)
    M2_true = [
        ε_i 0.0 0.0 0.0 0.0 0.0;
        0.0 ε_i 0.0 0.0 0.0 0.0;
        0.0 0.0 ε_i 0.0 0.0 0.0;
        0.0 0.0 0.0 μ_i 0.0 0.0;
        0.0 0.0 0.0 0.0 μ_i 0.0;
        0.0 0.0 0.0 0.0 0.0 μ_i
    ]

    M3 = TransferMatrix.construct_constitutive(ε)
    M3_true = [
        ε_i 0.0 0.0 0.0 0.0 0.0;
        0.0 ε_i 0.0 0.0 0.0 0.0;
        0.0 0.0 ε_i 0.0 0.0 0.0;
        0.0 0.0 0.0 1.0 0.0 0.0;
        0.0 0.0 0.0 0.0 1.0 0.0;
        0.0 0.0 0.0 0.0 0.0 1.0
    ]

    @test M1 == M1_true
    @test M2 == M2_true
    @test M3 == M3_true

    @testset "construct_constitutive full μ tensor" begin
        ε = SMatrix{3,3,ComplexF64}(2,0,0, 0,2,0, 0,0,2)
        μ = SMatrix{3,3,ComplexF64}(2,-0.6im,0, 0.6im,2,0, 0,0,1)
        M = TransferMatrix.construct_constitutive(ε, μ)
        @test M[4:6, 4:6] == μ
        @test M[1:3, 1:3] == ε
        @test all(==(0), M[1:3, 4:6]) && all(==(0), M[4:6, 1:3])

        ε_diag = Diagonal(SVector{3,ComplexF64}(2, 2, 2))
        M_bridge = TransferMatrix.construct_constitutive(ε_diag, μ)
        @test M_bridge == TransferMatrix.construct_constitutive(SMatrix{3,3,ComplexF64}(ε_diag), μ)
    end
end


@testset "construct_a" begin

    # Test orthorhombic crystal with principal axes
    # parallel to x, y, z. M is constant and diagonal.
    # The only nonzero coefficients are
    
    # a[3,5] = -ξ / M[3, 3]
    # a[6,2] =  ξ / M[6, 6]
    
    # Berreman, Optics in Stratefied and Anisotropic Media, 1972
    # DOI: 10.1364/JOSA.62.000502

    ξ = 15.0 + 15im

    ε = Diagonal([1 + 1im, 1 + 1im, 1 + 1im])
    μ = Diagonal([1 + 0im, 1 + 0im, 1 + 0im])
    ε[3,3] = 3.0 + 0im
    μ[3,3] = 0.0 + 5im

    M = @MMatrix zeros(ComplexF64, 6, 6)
    M[1:3, 1:3] = ε
    M[4:6, 4:6] = μ
    
    a = TransferMatrix.construct_a(ξ, M)

    to_subtract = @MMatrix zeros(ComplexF64, 6, 6)
    to_subtract[3,5] = -5.0 - 5im
    to_subtract[6,2] = 3.0 - 3im

    b = a - to_subtract
    test_against = zeros(ComplexF64, 6, 6)

    @test isapprox(b, zeros(ComplexF64, 6, 6), atol=1e-15)
end


@testset "construct_Δ" begin
 
    # Continue the test from `construct_a`.
    # For an orthorhombic crystal described above, the only nonzero elements of Δ are
    # Δ[2,1] = M[1,1] = ε[1,1]
    # Δ[4,3] = M[2,2] - ξ^2 / M[6,6] = ε[2,2] - ξ^2 / μ[3,3]
    # Δ[3,4] = M[4,4] = μ[1,1]
    # Δ[1,2] = M[5,5] - ξ^2 / M[3,3] = μ[2,2] - ξ^2 / ε[3,3]

    ξ = 0.

    ε = Diagonal([1 + 1im, 1 + 1im, 1 + 1im])
    μ = Diagonal([1 + 0im, 1 + 0im, 1 + 0im])
    ε[3,3] = 3.0 + 0im
    μ[3,3] = 0.0 + 5im

    M = zeros(ComplexF64, 6, 6)
    M[1:3, 1:3] = ε
    M[4:6, 4:6] = μ
    
    a = TransferMatrix.construct_a(ξ, M)
    Δ = TransferMatrix.construct_Δ(ξ, M, a)

    Δ21 = ε[1,1]
    Δ43 = ε[2,2] - ξ^2 / μ[3,3]
    Δ34 = μ[1,1]
    Δ12 = μ[2,2] - ξ^2 / ε[3,3]
    Δ_squared = Diagonal([Δ12 * Δ21, Δ12 * Δ21, Δ34 * Δ43, Δ34 * Δ43])
    Δ_cubed = Δ^3

    @test Δ^2 == Δ_squared
    @test Δ_cubed[1,1] == 0.0 + 0im
    @test Δ_cubed[1,2] == Δ12^2 * Δ21


    ξ = 1.0 + 1.0im
    M = rand(ComplexF64, 6, 6)
    a = rand(ComplexF64, 6, 6)
    Δ = TransferMatrix.construct_Δ(ξ, M, a)

    @test Δ[1,1] ==  M[5,1] + (M[5,3] + ξ) * a[3,1] + M[5,6] * a[6,1]
    @test Δ[1,2] ==  M[5,5] + (M[5,3] + ξ) * a[3,5] + M[5,6] * a[6,5]
    @test Δ[1,3] ==  M[5,2] + (M[5,3] + ξ) * a[3,2] + M[5,6] * a[6,2]
    @test Δ[1,4] == -M[5,4] - (M[5,3] + ξ) * a[3,4] - M[5,6] * a[6,4]
    @test Δ[2,1] ==  M[1,1] + M[1,3] * a[3,1] + M[1,6] * a[6,1]
    @test Δ[2,2] ==  M[1,5] + M[1,3] * a[3,5] + M[1,6] * a[6,5]
    @test Δ[2,3] ==  M[1,2] + M[1,3] * a[3,2] + M[1,6] * a[6,2]
    @test Δ[2,4] == -M[1,4] - M[1,3] * a[3,4] - M[1,6] * a[6,4]
    @test Δ[3,1] == -M[4,1] - M[4,3] * a[3,1] - M[4,6] * a[6,1]
    @test Δ[3,2] == -M[4,5] - M[4,3] * a[3,5] - M[4,6] * a[6,5]
    @test Δ[3,3] == -M[4,2] - M[4,3] * a[3,2] - M[4,6] * a[6,2]
    @test Δ[3,4] ==  M[4,4] + M[4,3] * a[3,4] + M[4,6] * a[6,4]
    @test Δ[4,1] ==  M[2,1] + M[2,3] * a[3,1] + (M[2,6] - ξ) * a[6,1]
    @test Δ[4,2] ==  M[2,5] + M[2,3] * a[3,5] + (M[2,6] - ξ) * a[6,5]
    @test Δ[4,3] ==  M[2,2] + M[2,3] * a[3,2] + (M[2,6] - ξ) * a[6,2]
    @test Δ[4,4] == -M[2,4] - M[2,3] * a[3,4] - (M[2,6] - ξ) * a[6,4]
end


@testset "calculate_E_modes" begin

    μ = 1.0

    ξ1 = 0.0
    ξ2 = 1.0
    ξ3 = 1.0 + 1im

    ε1 = I
    ε2 = fill(1, (3, 3))
    ε3 = [3 1 1 ; 1 2 1 ; 1 1 2]
    ε4 = [1 1 1 ; 1 1 1 ; 1 1 2]
    ε5 = [1 + 1im   1        1;
            1     1 + 1im    1;
            1       1      2 + 2im]

    q1 = [1., 1., 1., 1.]
    q2 = [1., 2., 1., 2.]
    q3 = [1. + 1im, 1. + 1im, 1. + 1im, 1. + 1im]
    q4 = [1. + 1im, 2. + 1im, 1. + 1im, 2. + 1im]

    γ1 = TransferMatrix.calculate_E_modes(ξ1, q1, ε1, μ)
    γ2 = TransferMatrix.calculate_E_modes(ξ1, q1, ε2, μ)
    γ3 = TransferMatrix.calculate_E_modes(ξ2, q2, ε3, μ)
    γ4 = TransferMatrix.calculate_E_modes(ξ2, q1, ε4, μ)
    γ5 = TransferMatrix.calculate_E_modes(ξ3, q3, ε5, μ)
    γ6 = TransferMatrix.calculate_E_modes(ξ3, q4, ε5, μ)

    γ1_test = ComplexF64[
        1 0 0;
        0 1 0;
       -1 0 0;
        0 1 0
        ]

    γ2_test = ComplexF64[
        1  0 -1;
        0  1 -1;
       -1  0  1;
        0  1 -1
        ]

    γ3_test = ComplexF64[
        1.0 -1.0 -1.0;
       -0.2  1.0 -0.4;
       -1.0  1.0  1.0;
       -0.2  1.0 -0.4
    ]

    γ4_test = ComplexF64[
        1.0  0.0 -2.0;
        0.0  1.0 -1.0;
       -1.0  0.0  2.0;
        0.0  1.0 -1.0
    ]

    γ5_test = ComplexF64[
        1.0  0.0 -0.5 - 1im;
        0.0  1.0 -0.5;
       -1.0  0.0  0.5 + 1im;
        0.0  1.0 -0.5
    ]

    γ6_test = ComplexF64[
          1.0+0.0im         -0.351351-0.108108im  -0.324324-0.945946im
    -0.166154+0.00923077im        1.0+0.0im           -0.32+0.24im
         -1.0+0.0im          0.351351+0.108108im   0.324324+0.945946im
    -0.166154+0.00923077im        1.0+0.0im           -0.32+0.24im
    ]

    for m in 1:4
        γ2_test[m, :] /= norm(γ2_test[m, :])
        γ3_test[m, :] /= norm(γ3_test[m, :])
        γ4_test[m, :] /= norm(γ4_test[m, :])
        γ5_test[m, :] /= norm(γ5_test[m, :])
        γ6_test[m, :] /= norm(γ6_test[m, :])
    end

    @test γ1 == γ1_test
    @test γ2 == γ2_test
    @test isapprox(γ3, γ3_test)
    @test isapprox(γ4, γ4_test)
    @test isapprox(γ5, γ5_test)
    @test isapprox(γ6, γ6_test, atol=1e-5)
end


@testset "dynamical_matrix" begin
    γ = ComplexF64[
         11 12 13;
         21 22 23;
         31 32 33;
         41 42 43]
    q = ComplexF64[1, 2, 3, 4]
    ξ = 1.0 + 1.0im
    μ = 2.0

    A_test = ComplexF64[
        11 21 31 41
        12 22 32 42
     -2 - 13im  19 - 23im    60 - 33im  121 - 43im
        12 44 96 168
    ]
    A_test[3, :] ./= μ
    A_test[4, :] ./= μ

    A = TransferMatrix.dynamical_matrix(ξ, q, γ, μ)

    for (i, col) in enumerate(eachrow(A))
        @test A[:, i] == A_test[:, i]
    end
end


@testset "propagation_matrix" begin
    q = [1., 2., 3., 4.]
    d = π / 2
    ω = c_0
    P_test = Diagonal([-1.0im, -1.0, 1.0im, 1.0])
    P = TransferMatrix.propagation_matrix(ω, q)

    @test P(d) ≈ P_test
end

@testset "Layer" begin
    # Isotropic layer construction
    iso_layer = Layer(λ -> 1.5, 0.1)
    @test iso_layer.thickness == 0.1
    @test !isanisotropic(iso_layer)

    # Anisotropic layer construction
    nx = λ -> 1.5
    ny = λ -> 1.6
    nz = λ -> 1.7
    aniso_layer = Layer(nx, ny, nz, 0.2)
    @test aniso_layer.thickness == 0.2
    @test isanisotropic(aniso_layer)
    @test aniso_layer.dispersion isa Tuple

    # get_refractive_indices for isotropic
    n_iso = get_refractive_indices(iso_layer, 1.0)
    @test n_iso == (1.5, 1.5, 1.5)

    # get_refractive_indices for anisotropic
    n_aniso = get_refractive_indices(aniso_layer, 1.0)
    @test n_aniso == (1.5, 1.6, 1.7)

    # Negative thickness should throw
    @test_throws DomainError Layer(λ -> 1.5, -0.1)
    @test_throws DomainError Layer(nx, ny, nz, -0.1)

    # Uniaxial construction (two indices equal)
    no = λ -> 1.658
    ne = λ -> 1.486
    uniaxial = Layer(no, no, ne, 0.5)
    @test isanisotropic(uniaxial)
    n_uni = get_refractive_indices(uniaxial, 1.0)
    @test n_uni[1] == n_uni[2]  # ordinary indices equal
    @test n_uni[1] != n_uni[3]  # extraordinary different

    # Rotated anisotropic layer
    rotated = Layer(no, no, ne, 0.5; euler=(0.0, π/4, 0.0))
    @test isanisotropic(rotated)
    @test isrotated(rotated)
    @test get_euler_angles(rotated) == (0.0, π/4, 0.0)

    # Unrotated layer should not be flagged as rotated
    @test !isrotated(uniaxial)
    @test get_euler_angles(uniaxial) == (0.0, 0.0, 0.0)

    # Isotropic layer should not be flagged as rotated
    @test !isrotated(iso_layer)
    @test get_euler_angles(iso_layer) == (0.0, 0.0, 0.0)
end

@testset "broadcastable Layer" begin
    layer = Layer(λ -> 1.5 + 0.0im, 0.1)
    λs = [1.0, 1.1, 1.2]

    # A single Layer is treated as a scalar under broadcasting
    @test Base.broadcastable(layer) isa Base.RefValue

    # Broadcasting a helper over one layer across many λ works without Ref
    # (previously threw MethodError: no method matching length(::Layer))
    bare = get_refractive_indices.(layer, λs)
    @test bare == get_refractive_indices.(Ref(layer), λs)
    @test length(bare) == length(λs)

    # A Vector{Layer} still broadcasts element-wise
    layers = [Layer(λ -> 1.0, 0.1), Layer(λ -> 2.0, 0.2)]
    @test getfield.(layers, :thickness) == [0.1, 0.2]
end

@testset "euler_rotation_matrix" begin
    # Identity rotation (no angles)
    R0 = euler_rotation_matrix(0.0, 0.0, 0.0)
    @test R0 ≈ I(3) atol=1e-14

    # Rotation matrix should be orthogonal: R * R' = I
    φ, θ, ψ = π/6, π/4, π/3
    R = euler_rotation_matrix(φ, θ, ψ)
    @test R * R' ≈ I(3) atol=1e-14
    @test R' * R ≈ I(3) atol=1e-14
    @test det(R) ≈ 1.0 atol=1e-14  # proper rotation (not reflection)

    # 90° rotation about z (φ=π/2, θ=0, ψ=0) should map x→y, y→-x
    Rz90 = euler_rotation_matrix(π/2, 0.0, 0.0)
    @test Rz90 * [1, 0, 0] ≈ [0, 1, 0] atol=1e-14
    @test Rz90 * [0, 1, 0] ≈ [-1, 0, 0] atol=1e-14
    @test Rz90 * [0, 0, 1] ≈ [0, 0, 1] atol=1e-14

    # 90° rotation about y (φ=0, θ=π/2, ψ=0) should map z→x, x→-z
    Ry90 = euler_rotation_matrix(0.0, π/2, 0.0)
    @test Ry90 * [0, 0, 1] ≈ [1, 0, 0] atol=1e-14
    @test Ry90 * [1, 0, 0] ≈ [0, 0, -1] atol=1e-14
    @test Ry90 * [0, 1, 0] ≈ [0, 1, 0] atol=1e-14
end

@testset "rotate_dielectric_tensor" begin
    # Rotation of isotropic tensor should give same tensor
    ε_iso = TransferMatrix.dielectric_tensor(2.0 + 0im, 2.0 + 0im, 2.0 + 0im)
    R = euler_rotation_matrix(π/6, π/4, π/3)
    ε_rot = rotate_dielectric_tensor(ε_iso, R)
    @test ε_rot ≈ SMatrix{3,3}(ε_iso) atol=1e-12

    # Identity rotation should preserve tensor
    ε_aniso = TransferMatrix.dielectric_tensor(1.5 + 0im, 2.0 + 0im, 2.5 + 0im)
    R_id = euler_rotation_matrix(0.0, 0.0, 0.0)
    ε_id = rotate_dielectric_tensor(ε_aniso, R_id)
    @test ε_id ≈ SMatrix{3,3}(ε_aniso) atol=1e-14

    # Rotated tensor should remain symmetric: ε_ij = ε_ji
    ε_rot2 = rotate_dielectric_tensor(ε_aniso, R)
    @test ε_rot2 ≈ ε_rot2' atol=1e-14

    # Trace should be preserved under rotation (invariant)
    @test tr(ε_rot2) ≈ tr(ε_aniso) atol=1e-12
end

@testset "per-layer μ" begin
    l0 = Layer(λ -> 1.5, 0.3)
    @test l0.mu === nothing
    @test TransferMatrix.ismagnetic(l0) == false
    @test TransferMatrix.get_permeability(l0, 1.0) === nothing

    ls = Layer(λ -> 1.5, 0.3; mu = 2.5)                 # isotropic magnetic
    @test TransferMatrix.ismagnetic(ls)
    @test TransferMatrix.get_permeability(ls, 1.0) == SMatrix{3,3,ComplexF64}(2.5*I)

    M = SMatrix{3,3,ComplexF64}(2,0,0, 0,2,0, 0,0,3)
    lm = Layer(λ -> 1.5, 0.3; mu = [2 0 0; 0 2 0; 0 0 3])  # constant tensor
    @test TransferMatrix.get_permeability(lm, 1.0) == M

    lf = Layer(λ -> 1.5, 0.3; mu = λ -> M)              # dispersive
    @test TransferMatrix.get_permeability(lf, 0.7) == M

    la = Layer(λ -> 1.6, λ -> 1.6, λ -> 1.8, 0.4; euler=(0.1,0.2,0.3), mu = 2.5)
    @test TransferMatrix.isanisotropic(la) && TransferMatrix.ismagnetic(la)
end

@testset "tensor-μ dynamical matrix reduces to scalar" begin
    ε = TransferMatrix.dielectric_tensor(4.0+0im, 6.25+0im, 9.0+0im)
    ξ = 0.5; μs = 1.3
    M = TransferMatrix.construct_constitutive(ε, TransferMatrix.permeability_tensor(μs,μs,μs))
    a = TransferMatrix.construct_a(ξ, M); Δ = TransferMatrix.construct_Δ(ξ, M, a)
    q, _ = TransferMatrix.calculate_q(Δ, a); q = ComplexF64.(q)
    γref = TransferMatrix.calculate_E_modes(ξ, q, ε, μs)
    Dref = TransferMatrix.dynamical_matrix(ξ, q, γref, μs)
    μmat = SMatrix{3,3,ComplexF64}(μs*I)
    γt = TransferMatrix.calculate_E_modes_tensor(ξ, q, ε, μmat)
    Dt = TransferMatrix.dynamical_matrix(ξ, q, γt, μmat)
    # D·P·D⁻¹ is basis-invariant ⇒ compare the layer transfer, not raw D
    P = TransferMatrix.propagation_matrix(2π, q)
    Tref = Matrix(Dref)*Matrix(P(0.4))*inv(Matrix(Dref))
    Tt   = Matrix(Dt)*Matrix(P(0.4))*inv(Matrix(Dt))
    @test maximum(abs, Tref .- Tt) < 1e-10
end
