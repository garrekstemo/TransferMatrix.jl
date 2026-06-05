using Test
using LinearAlgebra
using StaticArrays
using RefractiveIndex
using TransferMatrix

const Z0_test = sqrt(TransferMatrix.μ_0 / TransferMatrix.ε_0)

@testset "Sheet constructors" begin
    # Scalar constant -> diagonal isotropic tensor in SI Siemens
    s = TransferMatrix.Sheet(2.0e-4 + 1.0e-4im)
    σ = s.conductivity(1.0)
    @test σ isa SMatrix{2,2,ComplexF64}
    @test σ[1,1] == 2.0e-4 + 1.0e-4im
    @test σ[2,2] == 2.0e-4 + 1.0e-4im
    @test σ[1,2] == 0
    @test σ[2,1] == 0

    # Scalar function of λ
    sf = TransferMatrix.Sheet(λ -> 1.0e-4 / λ)
    @test sf.conductivity(2.0)[1,1] ≈ 5.0e-5
    @test sf.conductivity(2.0)[1,1] == sf.conductivity(2.0)[2,2]

    # Keyword anisotropic tensor (constants and functions mixed)
    st = TransferMatrix.Sheet(; xx = 1.0e-4, yy = 3.0e-4, xy = λ -> 1.0e-5)
    M = st.conductivity(1.0)
    @test M[1,1] ≈ 1.0e-4
    @test M[2,2] ≈ 3.0e-4
    @test M[1,2] ≈ 1.0e-5
    @test M[2,1] == 0

    # n,k + thickness -> σ = -i ω ε₀ d (n²-1), isotropic; agrees with closed form
    n = 4.0 + 0.2im
    d = 6.5e-4   # μm (~0.65 nm)
    λ = 0.6
    sheet_nk = TransferMatrix.Sheet(λ0 -> n, d)
    σ_expected = -im * (2π * TransferMatrix.c_0 / λ) * TransferMatrix.ε_0 * d * (n^2 - 1)
    @test sheet_nk.conductivity(λ)[1,1] ≈ σ_expected
    @test sheet_nk.conductivity(λ)[2,2] ≈ σ_expected
    @test sheet_nk.conductivity(λ)[1,2] == 0

    # In-plane anisotropic n,k
    sheet_aniso = TransferMatrix.Sheet(λ0 -> 4.0 + 0.0im, λ0 -> 5.0 + 0.0im, d)
    σx = -im * (2π * TransferMatrix.c_0 / λ) * TransferMatrix.ε_0 * d * ((4.0+0im)^2 - 1)
    σy = -im * (2π * TransferMatrix.c_0 / λ) * TransferMatrix.ε_0 * d * ((5.0+0im)^2 - 1)
    @test sheet_aniso.conductivity(λ)[1,1] ≈ σx
    @test sheet_aniso.conductivity(λ)[2,2] ≈ σy
end

@testset "sheet_matrix G structure" begin
    # G = G_phys^-1 in basis (Ex, Ey, H_y, -Hx): identity + (+σ̃) in rows 3,4.
    σxx = 1.0e-4 + 2.0e-4im
    σyy = 3.0e-4 - 1.0e-4im
    σxy = 0.5e-4im
    σyx = -0.2e-4im
    sheet = TransferMatrix.Sheet(; xx = σxx, yy = σyy, xy = σxy, yx = σyx)
    G = TransferMatrix.sheet_matrix(sheet, 1.0)

    @test G isa SMatrix{4,4,ComplexF64}
    # Top-left 2×2 identity, tangential E continuous
    @test G[1,1] == 1 && G[2,2] == 1
    @test G[1,2] == 0 && G[1,3] == 0 && G[1,4] == 0
    @test G[2,1] == 0 && G[2,3] == 0 && G[2,4] == 0
    # H rows carry +σ̃ = +Z₀ σ
    @test G[3,1] ≈ Z0_test * σxx
    @test G[3,2] ≈ Z0_test * σxy
    @test G[4,1] ≈ Z0_test * σyx
    @test G[4,2] ≈ Z0_test * σyy
    @test G[3,3] == 1 && G[4,4] == 1
    @test G[3,4] == 0 && G[4,3] == 0

    # Zero conductivity -> identity (no-op interface)
    G0 = TransferMatrix.sheet_matrix(TransferMatrix.Sheet(0.0 + 0.0im), 1.0)
    @test G0 ≈ SMatrix{4,4,ComplexF64}(I)
end

@testset "sheet σ=0 regression (transfer)" begin
    n1, n2, nf = 1.0, 1.5, 2.0
    d = 0.1
    air = Layer(λ -> complex(n1), 0.0)
    film = Layer(λ -> complex(nf), d)
    sub = Layer(λ -> complex(n2), 0.0)
    layers = [air, film, sub]

    base = transfer(0.6, layers; θ = 0.2)
    zero_sheet = Dict(1 => TransferMatrix.Sheet(0.0 + 0.0im))
    withσ0 = transfer(0.6, layers; θ = 0.2, sheets = zero_sheet)

    @test isapprox(base.Rpp, withσ0.Rpp; atol = 1e-12)
    @test isapprox(base.Rss, withσ0.Rss; atol = 1e-12)
    @test isapprox(base.Tpp, withσ0.Tpp; atol = 1e-12)
    @test isapprox(base.Tss, withσ0.Tss; atol = 1e-12)
end

@testset "analytic conductive interface (pins G sign)" begin
    # Single sheet between media n1 | n2; compare R to closed-form conductive Fresnel.
    n1, n2 = 1.0, 1.5
    σ_s = 2.0e-4 + 1.0e-4im
    g = Z0_test * σ_s                      # dimensionless σ̃
    air = Layer(λ -> complex(n1), 0.0)
    sub = Layer(λ -> complex(n2), 0.0)
    layers = [air, sub]
    sheets = Dict(1 => TransferMatrix.Sheet(σ_s))

    for θ in (0.0, π/6)
        c1 = cos(θ)
        c2 = sqrt(1 - (n1 / n2 * sin(θ))^2)        # real (n1 < n2)
        rs = (n1*c1 - n2*c2 - g) / (n1*c1 + n2*c2 + g)
        rp = (n2*c1 - n1*c2 + g*c1*c2) / (n2*c1 + n1*c2 + g*c1*c2)
        res = transfer(1.0, layers; θ = θ, sheets = sheets)
        @test isapprox(res.Rss, abs2(rs); atol = 1e-8)
        @test isapprox(res.Rpp, abs2(rp); atol = 1e-8)
    end
end

@testset "sweep_angle forwards sheets" begin
    n1, n2 = 1.0, 1.5
    σ_s = 1.5e-4 + 0.0im
    air = Layer(λ -> complex(n1), 0.0)
    sub = Layer(λ -> complex(n2), 0.0)
    layers = [air, sub]
    sheets = Dict(1 => TransferMatrix.Sheet(σ_s))
    λs = [1.0, 1.2]
    θs = [0.0, π/8]

    spec = sweep_angle(λs, θs, layers; sheets = sheets)
    for (ii, θ) in enumerate(θs), (jj, λ) in enumerate(λs)
        ref = transfer(λ, layers; θ = θ, sheets = sheets)
        @test isapprox(spec.Rpp[ii, jj], ref.Rpp; atol = 1e-12)
        @test isapprox(spec.Rss[ii, jj], ref.Rss; atol = 1e-12)
    end
    # And the sheet actually changes the result vs no sheet
    spec0 = sweep_angle(λs, θs, layers)
    @test !isapprox(spec.Rss[1, 1], spec0.Rss[1, 1]; atol = 1e-6)
end

@testset "_propagate_full and propagate contract" begin
    air = Layer(λ -> complex(1.0), 0.0)
    film = Layer(λ -> complex(1.5), 0.1)
    sub = Layer(λ -> complex(1.5), 0.0)
    layers = [air, film, sub]

    # propagate stays a 5-tuple (existing callers unaffected)
    out = TransferMatrix.propagate(0.6, layers)
    @test length(out) == 5
    Γ, S, Ds, Ps, γs = out
    @test length(Ds) == length(layers)

    # _propagate_full adds qs (6-tuple), one q-vector per layer
    full = TransferMatrix._propagate_full(0.6, layers)
    @test length(full) == 6
    qs = full[6]
    @test length(qs) == length(layers)
    @test length(qs[1]) == 4
end

@testset "efield refactor regression (no sheets)" begin
    air = Layer(λ -> complex(1.0), 0.0)
    film = Layer(λ -> complex(1.5), 0.1)
    sub = Layer(λ -> complex(1.0), 0.0)
    layers = [air, film, sub]
    ef = efield(0.6, layers; dz = 0.01)
    @test size(ef.p, 1) == 3
    @test size(ef.s, 1) == 3
    @test size(ef.p, 2) == length(ef.z)
    @test length(ef.boundaries) == 2
    # Tangential E continuous across the internal interface even with no sheet
    @test isapprox(ef.boundaries[1], 0.0; atol = 1e-12)
end

@testset "efield E continuity across a sheet" begin
    n0 = 1.0
    air = Layer(λ -> complex(n0), 0.0)
    spacerL = Layer(λ -> complex(1.0), 0.5)
    spacerR = Layer(λ -> complex(1.0), 0.5)
    sub = Layer(λ -> complex(n0), 0.0)
    layers = [air, spacerL, spacerR, sub]
    sheets = Dict(2 => TransferMatrix.Sheet(2.0e-4 + 1.0e-4im))   # interface between layers 2 and 3
    # A sheet introduces a jump in tangential H (surface current) but NOT in
    # tangential E, so the straddling E samples differ only by the field's own
    # propagation across the one-dz gap (an O(dz) effect that vanishes as dz→0,
    # unlike a real discontinuity). dz is chosen small enough that this gap stays
    # well under the continuity tolerance.
    dz = 1.0e-4
    ef = efield(0.6, layers; dz = dz, sheets = sheets)

    # Find samples straddling the layer-2/3 interface (at z = boundaries[2])
    zint = ef.boundaries[2]
    jbelow = findlast(z -> z ≤ zint, ef.z)
    jabove = jbelow + 1
    # Tangential E (Ex, Ey) continuous across the sheet
    @test isapprox(ef.p[1, jbelow], ef.p[1, jabove]; atol = 5e-3)   # Ex (p)
    @test isapprox(ef.s[2, jbelow], ef.s[2, jabove]; atol = 5e-3)   # Ey (s)
end

@testset "hfield basics and shared z-grid" begin
    air = Layer(λ -> complex(1.0), 0.0)
    film = Layer(λ -> complex(2.0), 0.2)
    sub = Layer(λ -> complex(1.5), 0.0)
    layers = [air, film, sub]

    hf = hfield(0.6, layers; dz = 0.01)
    ef = efield(0.6, layers; dz = 0.01)
    @test hf isa MagneticField
    @test size(hf.p, 1) == 3
    @test size(hf.s, 1) == 3
    @test size(hf.p, 2) == length(hf.z)
    @test hf.z == ef.z                       # shared grid
    @test hf.boundaries == ef.boundaries
end

@testset "transfer vs field cross-check (with sheet)" begin
    # Guard against a placement/sign mismatch between the two sheet-injection
    # code paths: transfer's forward _propagate_core recursion vs _field's
    # backward mode-coefficient recursion (used by efield/hfield). They are
    # maintained separately, so we reconstruct R from the FIELD path and require
    # it to match transfer's R.
    air = Layer(λ -> complex(1.0), 0.0)
    spacer = Layer(λ -> complex(1.0), 0.5)
    sub = Layer(λ -> complex(1.5), 0.0)
    layers = [air, spacer, sub]
    sheets = Dict(2 => TransferMatrix.Sheet(2.0e-4 + 1.0e-4im))

    res = transfer(0.6, layers; θ = 0.1, sheets = sheets)

    # _field normalizes the substrate to a transmitted-only drive, so its layer-1
    # mode coefficients (sampled at the incident-side z) reproduce the incident and
    # reflected amplitudes. Mode order is (p-trans, s-trans, p-refl, s-refl):
    # the forward modes carry the unit incident drive, the backward modes carry r.
    F = TransferMatrix._field(0.6, layers; θ = 0.1, sheets = sheets)
    @test F.layer_of_z[1] == 1                       # first sample is inside the incident medium
    ap = F.amp_p[:, 1]
    as = F.amp_s[:, 1]
    rpp_field = ap[3] / ap[1]                         # p-refl / p-incident
    rss_field = as[4] / as[2]                         # s-refl / s-incident
    @test isapprox(res.Rpp, abs2(rpp_field); atol = 1e-10)
    @test isapprox(res.Rss, abs2(rss_field); atol = 1e-10)
end

@testset "hfield H-jump across a sheet" begin
    # In-plane H jumps by σ̃·E across the sheet while tangential E stays continuous.
    n0 = 1.0
    air = Layer(λ -> complex(n0), 0.0)
    spacerL = Layer(λ -> complex(1.0), 0.5)
    spacerR = Layer(λ -> complex(1.0), 0.5)
    sub = Layer(λ -> complex(n0), 0.0)
    layers = [air, spacerL, spacerR, sub]
    σ_s = 3.0e-4 + 0.0im
    g = Z0_test * σ_s
    sheets = Dict(2 => TransferMatrix.Sheet(σ_s))
    dz = 5e-4
    ef = efield(0.6, layers; dz = dz, sheets = sheets)
    hf = hfield(0.6, layers; dz = dz, sheets = sheets)

    zint = ef.boundaries[2]
    jb = findlast(z -> z ≤ zint, ef.z)
    ja = jb + 1
    # s-incidence: E is along y; expected ΔHx = σ̃ Ey (use field at the interface)
    Ey = ef.s[2, jb]
    ΔHx = hf.s[1, ja] - hf.s[1, jb]
    @test isapprox(ΔHx, g * Ey; rtol = 0.05, atol = 1e-6)
end
