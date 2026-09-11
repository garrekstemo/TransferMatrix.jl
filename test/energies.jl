# mode_amplitudes, layer_energies and amplitude_coefficients (v3.4.0): the
# per-layer Berreman amplitudes behind efield, the closed-form ∫ E†εE dz over
# each finite layer, and the public r/t amplitude coefficients.

# Midpoint-rule ∫ E†εH E dz over each finite layer from an efield profile. The
# stacks below give the incident layer a thickness of (K + 1/2)·dz and every
# other layer an integer multiple of dz, so all samples sit half a step away
# from every interface and each layer is integrated by exactly d/dz interior
# samples (error O(dz²), no interface sample ambiguity for discontinuous Ez).
function midpoint_layer_energies(ef, layers, λ, dz)
    N = length(layers)
    bnds = ef.boundaries
    Up = zeros(N - 2)
    Us = zeros(N - 2)
    for i in 2:N-1
        ε = SMatrix{3,3,ComplexF64}(TransferMatrix._layer_epsilon(layers[i], λ))
        εH = (ε + ε') / 2
        idx = findall(z -> bnds[i-1] < z < bnds[i], ef.z)
        @test length(idx) == round(Int, layers[i].thickness / dz)
        for j in idx
            Ep = SVector{3,ComplexF64}(ef.p[1, j], ef.p[2, j], ef.p[3, j])
            Es = SVector{3,ComplexF64}(ef.s[1, j], ef.s[2, j], ef.s[3, j])
            Up[i-1] += real(dot(Ep, εH * Ep)) * dz
            Us[i-1] += real(dot(Es, εH * Es)) * dz
        end
    end
    return Up, Us
end

# Field at z inside layer i reconstructed from the entry-face amplitudes.
function modal_field(A, i, z)
    zref = i == 1 ? 0.0 : A.boundaries[i-1]
    k0 = 2π
    E = zeros(ComplexF64, 3)
    for m in 1:4
        E .+= A.p[i][m] * exp(im * k0 * A.qs[i][m] * (z - zref)) * A.E_modes_per_layer[i][m, :]
    end
    return E
end

@testset "mode amplitudes and layer energies" begin

    @testset "isotropic DBR cavity, normal incidence" begin
        λ = 1.0
        dz = 5e-5
        nH = 2.3; nL = 1.45; ns = 1.5
        dH = 0.1087; dL = 0.1724; dsp = 0.3334              # integer multiples of dz
        pair = [Layer(λ -> nH, dH), Layer(λ -> nL, dL)]
        mirror = vcat([pair for _ in 1:4]...)
        layers = [Layer(λ -> 1.0, 0.1 + dz / 2); mirror; Layer(λ -> ns, dsp); reverse(mirror); Layer(λ -> 1.45, 0.3)]
        N = length(layers)

        A = mode_amplitudes(λ, layers)
        c = amplitude_coefficients(λ, layers; method=:eig)
        @test length(A.p) == N && length(A.s) == N
        @test length(A.boundaries) == N - 1 && A.boundaries[1] == 0.0
        # Incident medium referenced to z = 0: (1, 0, rpp, rps); exit medium at
        # the last interface: (tpp, tps, 0, 0). Only slots 1/3 (p) and 2/4 (s)
        # are populated for an isotropic layer at normal incidence.
        @test A.p[1] ≈ SVector(1.0, 0.0, c.rpp, c.rps) atol=1e-12
        @test A.s[1] ≈ SVector(0.0, 1.0, c.rsp, c.rss) atol=1e-12
        @test A.p[end] ≈ SVector(c.tpp, c.tps, 0.0, 0.0) atol=1e-12
        @test A.s[end] ≈ SVector(c.tsp, c.tss, 0.0, 0.0) atol=1e-12
        for i in 1:N
            @test A.p[i][2] == 0 && A.p[i][4] == 0
            @test A.s[i][1] == 0 && A.s[i][3] == 0
            # forward slots have Re q > 0, backward slots Re q < 0
            @test real(A.qs[i][1]) > 0 && real(A.qs[i][2]) > 0
            @test real(A.qs[i][3]) < 0 && real(A.qs[i][4]) < 0
        end

        # Reference-plane check: the modal expansion from the entry face
        # reproduces efield at arbitrary interior points of every layer.
        ef = efield(λ, layers; dz=dz)
        for i in 1:N
            zl = i == 1 ? -layers[1].thickness : A.boundaries[i-1]
            zr = i == N ? A.boundaries[end] + layers[end].thickness : A.boundaries[i]
            for frac in (0.13, 0.5, 0.87)
                z = zl + frac * (zr - zl)
                j = argmin(abs.(ef.z .- z))
                @test modal_field(A, i, ef.z[j]) ≈ ef.p[:, j] atol=1e-10
            end
        end

        # Closed form vs midpoint-rule integration of the sampled field.
        U = layer_energies(λ, layers)
        Up, Us = midpoint_layer_energies(ef, layers, λ, dz)
        @test length(U.p) == N - 2 && length(U.s) == N - 2
        for j in 1:N-2
            @test isapprox(U.p[j], Up[j]; rtol=1e-6)
            @test isapprox(U.s[j], Us[j]; rtol=1e-6)
            @test U.p[j] > 0 && U.s[j] > 0
        end
        # p and s are physically identical at normal incidence.
        @test U.p ≈ U.s rtol=1e-12
        # Units: a unit-index layer returns ∫|E|² dz — check with a bare slab.
        slab = [Layer(λ -> 1.0, 0.1 + dz / 2), Layer(λ -> 1.0, 0.3), Layer(λ -> 1.0, 0.2)]
        @test layer_energies(λ, slab).p[1] ≈ 0.3 rtol=1e-12
    end

    @testset "regression vs. scalar 2×2 solver (CaF2 | DBR | gap | DBR | CaF2)" begin
        # Cavity from the cavity-feedback-loop analysis: CaF2 substrate,
        # 124.5 nm SiO2 + 10 alternating Ge/ZnS layers (ZnS facing the gap),
        # gap n = 1.43 of length L, the mirror reversed, CaF2. Ge and ZnS indices
        # carry the calibration offsets +0.037 / −0.07; the |λ| clamps keep the
        # Sellmeier fits inside their valid ranges (inactive at the test λ).
        mirror_spec = [
            (:SiO2, 0.1245),
            (:Ge, 0.22356), (:ZnS, 0.58645),
            (:Ge, 0.27087), (:ZnS, 0.58645),
            (:Ge, 0.27087), (:ZnS, 0.58645),
            (:Ge, 0.27087), (:ZnS, 0.99514),
            (:Ge, 0.22383), (:ZnS, 0.51682),
        ]
        n_ge = RefractiveMaterial("main", "Ge", "Icenogle")
        n_zns = RefractiveMaterial("main", "ZnS", "Amotchkina")[1]   # [1] is the n fit
        n_sio2 = RefractiveMaterial("main", "SiO2", "Malitson")
        n_caf2 = RefractiveMaterial("main", "CaF2", "Malitson")
        mirror_index(mat) =
            mat == :Ge ? (λ -> real(n_ge(clamp(abs(λ), 2.5, 12.0))) + 0.037) :
            mat == :ZnS ? (λ -> real(n_zns(clamp(abs(λ), 0.41, 13.9))) - 0.07) :
            mat == :SiO2 ? (λ -> real(n_sio2(clamp(abs(λ), 0.22, 6.6)))) :
            (λ -> real(n_caf2(clamp(abs(λ), 0.24, 9.6))))
        mirror = [Layer(mirror_index(mat), t) for (mat, t) in mirror_spec]
        cavity(L, n) = [Layer(mirror_index(:CaF2), 1.0); mirror; Layer(λ -> n, L);
                        reverse(mirror); Layer(mirror_index(:CaF2), 1.0)]
        gap = length(mirror) + 2
        n_gap = 1.43

        # Guard against a drifted refractive-index database: the indices the
        # reference solver saw at the m = 3 resonance (printed to 8 decimals).
        λ3 = 1e4 / 1896.2214732815
        @test mirror_index(:CaF2)(λ3) ≈ 1.39557766 atol=1e-7
        @test mirror_index(:SiO2)(λ3) ≈ 1.32229921 atol=1e-7
        @test mirror_index(:Ge)(λ3) ≈ 4.05162911 atol=1e-7
        @test mirror_index(:ZnS)(λ3) ≈ 2.18552332 atol=1e-7

        # (m, L_um, ν_cm⁻¹, f_E) from the reference solver: f_E = U_gap / ΣU with
        # U = ∫ n²|E|² dz over each finite layer, at the m-th gap resonance
        # nearest 1771.3 cm⁻¹ (L = m / (2 n ν₀)).
        baseline = [
            (1,  1.9739758914, 1949.9623300621, 0.224916651464),
            (2,  3.9479517829, 1918.8390492418, 0.384452915412),
            (3,  5.9219276743, 1896.2214732815, 0.496381620887),
            (4,  7.8959035658, 1879.2842336809, 0.576742184625),
            (5,  9.8698794572, 1866.2359922239, 0.636263649840),
            (6, 11.8438553486, 1855.9109699046, 0.681736121937),
            (7, 13.8178312401, 1847.5689428509, 0.717407261587),
            (8, 15.7918071315, 1840.7041923499, 0.746034941529),
            (9, 17.7657830229, 1834.9522608400, 0.769497800740),
        ]
        for (m, L, ν, f_ref) in baseline
            U = layer_energies(1e4 / ν, cavity(L, n_gap))
            @test length(U.p) == 2 * length(mirror) + 1
            @test isapprox(U.p[gap-1] / sum(U.p), f_ref; rtol=1e-8)
            @test isapprox(U.s[gap-1] / sum(U.s), f_ref; rtol=1e-8)
        end

        # Per-layer amplitudes at m = 3. The reference solver writes
        # E = a e^{ikz} + b e^{−ikz} with z from the layer's LEFT face (layer 1
        # from the first interface) and unit transmitted amplitude a[end] = 1.
        # Same reference planes as ModeAmplitudes, so the mapping is
        #   s incidence: (a, b) = c · (A.s[2],  A.s[4])
        #   p incidence: (a, b) = c · (A.p[1], −A.p[3])   (backward p vector is −x̂)
        # with one overall scale c = a[1] (the package uses unit incident
        # amplitude, the reference unit transmitted amplitude).
        amps = [   # (Re a, Im a, Re b, Im b)
            (-9.976884934903e-01, -6.799520737226e-02,  5.551115123126e-16, -2.382911781413e-03),
            (-1.025333171198e+00, -6.981324053305e-02,  2.764467770816e-02, -5.648786206252e-04),
            (-6.488250634392e-01, -1.799309820372e-01, -3.162382166188e-01, -9.430794757296e-02),
            (-1.114719120453e-01, -1.037391886343e+00, -2.686073882455e-01,  6.147378388986e-01),
            ( 9.327031075782e-01, -5.264412730738e-02,  7.012275742150e-01,  1.911487838531e-01),
            ( 2.621528484831e-01,  1.533325155537e+00,  3.984115137754e-01, -1.273782964571e+00),
            (-1.459311001850e+00,  1.487042627038e-01, -1.316229259116e+00, -2.734841360962e-01),
            (-4.877718118763e-01, -2.467568984780e+00, -6.421040556675e-01,  2.296946896071e+00),
            ( 2.403072048665e+00, -2.872724391585e-01,  2.307483287906e+00,  4.339966345394e-01),
            ( 8.528233902418e-01,  4.107182800141e+00,  1.069346047090e+00, -3.976707425354e+00),
            (-2.901671849542e+00, -1.699458367329e+00, -2.964414380945e+00,  1.475010888142e+00),
            ( 2.299711493285e-01, -6.207608229247e+00, -1.918716181529e-01,  6.157265348030e+00),
            ( 6.140223054896e+00, -1.880180121103e+00,  5.921351071507e+00,  2.280151202158e+00),
            (-6.005618313063e+00, -1.587528109423e+00, -6.054318427028e+00,  1.137703753671e+00),
            ( 1.038393502038e-01, -3.361130491230e+00, -1.385477020618e-01,  3.308222604566e+00),
            ( 3.082690571232e+00, -2.845049054554e+00,  2.784562677095e+00,  3.033964502870e+00),
            (-1.058325778598e+00,  2.176539011198e+00, -8.679613404246e-01, -2.181643900783e+00),
            (-2.395993495063e+00, -7.653973798387e-01, -2.313300273460e+00,  5.802623965298e-01),
            ( 6.181391797893e-01, -1.330149098977e+00,  5.203198792502e-01,  1.239439942772e+00),
            ( 1.493433670101e+00,  4.345569072845e-01,  1.283878005420e+00, -3.637110479852e-01),
            (-3.546515212566e-01,  8.639906036962e-01, -3.214151011304e-01, -6.515376235892e-01),
            (-1.017781175247e+00, -2.279372372893e-01, -6.214602736051e-01,  2.511287712077e-01),
            ( 1.925432122656e-01, -6.448092117488e-01,  2.143888345420e-01,  2.498813981698e-01),
            ( 1.008003414529e+00, -2.002856536256e-01, -2.717743881319e-02, -5.400032398810e-03),
            ( 1.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00),
        ]
        a_ref = [ComplexF64(x[1], x[2]) for x in amps]
        b_ref = [ComplexF64(x[3], x[4]) for x in amps]
        L3 = 3 / (2 * n_gap * 1771.3e-4)
        A = mode_amplitudes(λ3, cavity(L3, n_gap))
        @test length(A.p) == length(amps)
        c_s = a_ref[1] / A.s[1][2]
        c_p = a_ref[1] / A.p[1][1]
        @test abs(c_s) ≈ 1 / abs(amplitude_coefficients(λ3, cavity(L3, n_gap)).tss) rtol=1e-8
        for i in eachindex(amps)
            @test isapprox(c_s * A.s[i][2], a_ref[i]; atol=1e-8)
            @test isapprox(c_s * A.s[i][4], b_ref[i]; atol=1e-8)
            @test isapprox(c_p * A.p[i][1], a_ref[i]; atol=1e-8)
            @test isapprox(-c_p * A.p[i][3], b_ref[i]; atol=1e-8)
        end
    end

    @testset "anisotropic layer at oblique incidence" begin
        λ = 1.0
        dz = 5e-5
        θ = 0.6
        uni = Layer(λ -> 1.6, λ -> 1.6, λ -> 1.9, 0.25; euler=(0.3, 0.7, 0.2))
        layers = [Layer(λ -> 1.0, 0.1 + dz / 2), Layer(λ -> 1.5, 0.2), uni,
                  Layer(λ -> 1.5 + 0.05im, 0.15), Layer(λ -> 1.45, 0.3)]
        N = length(layers)
        ef = efield(λ, layers; θ=θ, dz=dz)
        U = layer_energies(λ, layers; θ=θ)
        Up, Us = midpoint_layer_energies(ef, layers, λ, dz)
        for j in 1:N-2
            @test isapprox(U.p[j], Up[j]; rtol=1e-5)
            @test isapprox(U.s[j], Us[j]; rtol=1e-5)
        end
        # The rotated crystal mixes polarizations: all four slots are populated.
        A = mode_amplitudes(λ, layers; θ=θ)
        @test all(abs.(A.p[3]) .> 1e-3)
        @test all(abs.(A.s[3]) .> 1e-3)
        @test A.k_par ≈ sin(θ)
    end

    @testset "evanescent layer (degenerate phase integral)" begin
        # Frustrated TIR: glass | air gap | glass beyond the critical angle. The
        # gap modes are q = ±iκ, so the forward/backward cross terms have
        # q_m' − conj(q_m) = 0 and the integrand is constant in z.
        λ = 1.0
        dz = 5e-5
        θ = 1.1
        layers = [Layer(λ -> 1.5, 0.1 + dz / 2), Layer(λ -> 1.0, 0.3), Layer(λ -> 1.5, 0.3)]
        A = mode_amplitudes(λ, layers; θ=θ)
        @test abs(real(A.qs[2][1])) < 1e-12 && imag(A.qs[2][1]) > 0
        @test abs(real(A.qs[2][3])) < 1e-12 && imag(A.qs[2][3]) < 0
        ef = efield(λ, layers; θ=θ, dz=dz)
        U = layer_energies(λ, layers; θ=θ)
        Up, Us = midpoint_layer_energies(ef, layers, λ, dz)
        @test isapprox(U.p[1], Up[1]; rtol=1e-6)
        @test isapprox(U.s[1], Us[1]; rtol=1e-6)
    end

    @testset "Unitful and sheet arguments pass through" begin
        layers = [Layer(λ -> 1.0, 0.1), Layer(λ -> 1.5, 0.3), Layer(λ -> 1.45, 0.3)]
        @test mode_amplitudes(1u"μm", layers).p == mode_amplitudes(1.0, layers).p
        @test layer_energies(1000u"nm", layers).p == layer_energies(1.0, layers).p
        sheets = Dict(1 => Sheet(6.08e-5 + 1.0e-5im))
        U0 = layer_energies(1.0, layers)
        U1 = layer_energies(1.0, layers; sheets=sheets)
        @test !(U0.p ≈ U1.p)
        # sheet-modified amplitudes still reproduce efield
        ef = efield(1.0, layers; dz=1e-3, sheets=sheets)
        A = mode_amplitudes(1.0, layers; sheets=sheets)
        j = argmin(abs.(ef.z .- 0.15))
        @test modal_field(A, 2, ef.z[j]) ≈ ef.p[:, j] atol=1e-10
        @test_throws ArgumentError mode_amplitudes(1.0, layers; sheets=Dict(3 => Sheet(1e-4 + 0im)))
    end
end

@testset "amplitude_coefficients" begin

    @testset "single interface Fresnel amplitudes" begin
        layers = [Layer(λ -> 1.0, 0.1), Layer(λ -> 1.5, 0.1)]
        c = amplitude_coefficients(1.0, layers)
        # rss is the tangential-E reflection coefficient; rpp carries the −x̂
        # sign of the backward p mode vector, so rpp = −rss at normal incidence.
        @test c.rss ≈ (1.0 - 1.5) / (1.0 + 1.5) atol=1e-12
        @test c.rpp ≈ -c.rss atol=1e-12
        @test c.tss ≈ 2 / 2.5 atol=1e-12
        @test c.tpp ≈ c.tss atol=1e-12
        @test c.rps == 0 && c.rsp == 0 && c.tps == 0 && c.tsp == 0
        @test keys(c) == (:rpp, :rps, :rsp, :rss, :tpp, :tps, :tsp, :tss)
    end

    @testset "|r|², |t|² vs transfer for a lossless stack" begin
        λ = 1.0
        air = Layer(λ -> 1.0, 0.1)
        films = [Layer(λ -> 1.5, 0.3), Layer(λ -> 2.0, 0.2), Layer(λ -> 1.7, 0.15)]
        for θ in (0.0, 0.4, 0.9)
            # same medium on both sides: T = |t|²
            layers = [air; films; air]
            c = amplitude_coefficients(λ, layers; θ=θ)
            R = transfer(λ, layers; θ=θ)
            @test abs2(c.rpp) ≈ R.Rpp atol=1e-12
            @test abs2(c.rss) ≈ R.Rss atol=1e-12
            @test abs2(c.tpp) ≈ R.Tpp atol=1e-12
            @test abs2(c.tss) ≈ R.Tss atol=1e-12
            @test abs2(c.rpp) + abs2(c.tpp) ≈ 1 atol=1e-10
            # different exit medium: T = |t|² Re(q_out)/q_in with unit-norm modes
            n_out = 1.45
            layers2 = [air; films; Layer(λ -> n_out, 0.1)]
            c2 = amplitude_coefficients(λ, layers2; θ=θ)
            R2 = transfer(λ, layers2; θ=θ)
            q_in = cos(θ)
            q_out = sqrt(n_out^2 - sin(θ)^2)
            @test abs2(c2.rpp) ≈ R2.Rpp atol=1e-12
            @test abs2(c2.rss) ≈ R2.Rss atol=1e-12
            @test abs2(c2.tpp) * q_out / q_in ≈ R2.Tpp atol=1e-12
            @test abs2(c2.tss) * q_out / q_in ≈ R2.Tss atol=1e-12
            # both backends agree
            ce = amplitude_coefficients(λ, layers2; θ=θ, method=:eig)
            for k in keys(c2)
                @test isapprox(c2[k], ce[k]; atol=1e-12)
            end
        end
        @test_throws ArgumentError amplitude_coefficients(λ, [air, air]; method=:bogus)
    end

    @testset "perfect-conductor limit: rss → −1, tss → 0" begin
        # Thick gold-like layer (n ≈ 0.9 + 20i near 3 μm) behind a CaF2 medium.
        λ = 1e4 / 3000
        gold = Layer(λ -> 0.9 + 20.0im, 2.0)
        layers = [Layer(λ -> 1.40, 1.0), gold, Layer(λ -> 1.0, 1.0)]
        c = amplitude_coefficients(λ, layers)
        @test abs(c.rss + 1) < 0.15
        @test real(c.rss) < -0.95
        @test abs(c.tss) < 1e-10
        @test c.rpp ≈ -c.rss atol=1e-12
        @test c.tpp ≈ c.tss atol=1e-12
        @test abs2(c.rss) ≈ transfer(λ, layers).Rss atol=1e-12
        # the bulk-metal Fresnel value, for reference
        n1 = 1.40; n2 = 0.9 + 20.0im
        @test c.rss ≈ (n1 - n2) / (n1 + n2) atol=1e-6
    end
end
