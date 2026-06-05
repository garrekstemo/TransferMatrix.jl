# 2D Conductive Sheets (TMDC monolayers) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add zero-thickness anisotropic conductive sheets to the Berreman 4×4 TMM so a TMDC monolayer can be simulated inside a cavity via its surface conductivity, with R/T spectra, sweeps, and E/H field profiles.

**Architecture:** A new `Sheet` type holds a 2×2 SI sheet-conductivity tensor `σ(λ)`. `sheet_matrix` builds the 4×4 interface matrix `G = G_phys⁻¹` (`+σ̃`, `σ̃ = Z₀σ`) in the dynamical-matrix basis `(Eₓ,Eᵧ,H_y,−Hₓ)`. The transfer product is reformulated to interface–propagation (L–P) form so a sheet injects at one interface as `Dᵢ⁻¹ G Dᵢ₊₁`. Sheets are passed as a side-channel keyword `sheets=Dict(i => sheet)` (sheet at interface between layers `i` and `i+1`). `efield`/`hfield` share a `_field` core; H reconstruction reuses `dynamical_matrix` rows 3–4 plus `H_z`.

**Tech Stack:** Julia 1.10 LTS, StaticArrays, RefractiveIndex, LinearAlgebra, Test + Aqua.

**Spec:** `docs/superpowers/specs/2026-06-05-conductive-sheets-design.md`

**Conventions (from CLAUDE.md):** single `#` comments only (never `##`); module-level `const` is fine (only avoid `const` in throwaway scripts); prefer `eachindex`; run tests from project root with `julia --project=. -e 'using Pkg; Pkg.test()'`.

**Deviation from spec §4:** `sheet_matrix` lives in `src/sheet.jl` (with the `Sheet` type), not `matrix_constructors.jl`. Reason: `sheet_matrix(sheet::Sheet, λ)`'s signature needs `Sheet` defined first, and co-locating all sheet code in one file is cleaner. `sheet.jl` is `include`d after `layer.jl` (so `refractive_index` exists) and before `general_TMM.jl`.

**Sheet default sentinel:** Public functions default `sheets=nothing` (a true no-sheets sentinel) rather than an empty `Dict`, so the no-sheet hot path allocates nothing. Inputs (`Dict` or iterable of `i => sheet` pairs) are normalized to `Dict{Int,Sheet}` at the API boundary.

---

## Task 1: `Sheet` type, constructors, module wiring

**Files:**
- Create: `src/sheet.jl`
- Modify: `src/TransferMatrix.jl` (add `include("sheet.jl")`; export `Sheet`)
- Test: `test/sheets.jl` (new), `test/runtests.jl` (register)

- [ ] **Step 1: Create the test file and register it**

Create `test/sheets.jl`:

```julia
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
```

Register in `test/runtests.jl` by adding after line 18 (`include("integration.jl")`):

```julia
include("sheets.jl")
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `UndefVarError: Sheet not defined` (type does not exist yet).

- [ ] **Step 3: Create `src/sheet.jl`**

```julia
"""
    Sheet(σ)
    Sheet(; xx, yy, xy=0, yx=0)
    Sheet(material, d)
    Sheet(nx, ny, d)

A zero-thickness 2D conductive sheet (e.g. a TMDC monolayer) for the transfer
matrix. The internal representation is a callable `λ -> SMatrix{2,2,ComplexF64}`
returning the **SI sheet conductivity tensor** `[σxx σxy; σyx σyy]` in Siemens.

- `Sheet(σ)` — scalar isotropic conductivity (`σ::Number` or `σ(λ)->Number`),
  stored as `diag(σ, σ)`.
- `Sheet(; xx, yy, xy=0, yx=0)` — anisotropic tensor; each entry a `Number` or
  `λ->Number`.
- `Sheet(material, d)` — convert a refractive index to a sheet via
  `σ = -i ω ε₀ d (n²-1)` (isotropic). `material` is a `RefractiveMaterial` or
  `n(λ)`; `d` is the effective thickness in μm.
- `Sheet(nx, ny, d)` — in-plane anisotropic index conversion (`σxx` from `nx`,
  `σyy` from `ny`).

The factor `Z₀ = √(μ₀/ε₀)` is applied exactly once, later, in [`sheet_matrix`](@ref).
"""
struct Sheet{F}
    conductivity::F
    Sheet{F}(c::F) where {F} = new{F}(c)
end

# Internal: store a raw λ -> SMatrix closure without scalar-wrapping.
_rawsheet(f) = Sheet{typeof(f)}(f)

# Internal: evaluate a component that may be a constant or a function of λ.
_cval(x::Number, λ) = ComplexF64(x)
_cval(f, λ) = ComplexF64(f(λ))

# Internal: turn an index spec into an n(λ) function.
_index_fn(m::RefractiveMaterial) = refractive_index(m)
_index_fn(x::Number) = _ -> x
_index_fn(f::Function) = f

# Internal: SI scalar sheet conductivity from an index function and thickness d (μm).
_sigma_from_index(nfun, d) = λ -> -im * (2π * c_0 / λ) * ε_0 * d * (ComplexF64(nfun(λ))^2 - 1)

# Internal: wrap a scalar conductivity (constant or function) as a diagonal sheet.
_diagonal_sheet(σ) = _rawsheet(λ -> (s = _cval(σ, λ); SMatrix{2,2,ComplexF64}(s, 0, 0, s)))

# Scalar isotropic conductivity (Number or function returning a Number).
Sheet(σ::Number) = _diagonal_sheet(σ)
Sheet(σ::Function) = _diagonal_sheet(σ)

# Anisotropic tensor by component. SMatrix is column-major: (xx, yx, xy, yy).
function Sheet(; xx, yy, xy = 0, yx = 0)
    _rawsheet(λ -> SMatrix{2,2,ComplexF64}(_cval(xx, λ), _cval(yx, λ), _cval(xy, λ), _cval(yy, λ)))
end

# Index -> conductivity (isotropic).
Sheet(material::RefractiveMaterial, d::Real) = _diagonal_sheet(_sigma_from_index(_index_fn(material), d))
Sheet(n::Function, d::Real) = _diagonal_sheet(_sigma_from_index(n, d))

# Index -> conductivity (in-plane anisotropic).
function Sheet(nx, ny, d::Real)
    σx = _sigma_from_index(_index_fn(nx), d)
    σy = _sigma_from_index(_index_fn(ny), d)
    _rawsheet(λ -> SMatrix{2,2,ComplexF64}(σx(λ), 0, 0, σy(λ)))
end
```

- [ ] **Step 4: Wire `sheet.jl` into the module and export `Sheet`**

In `src/TransferMatrix.jl`, change the include block (currently lines 36–39) to:

```julia
include("matrix_constructors.jl")
include("layer.jl")
include("sheet.jl")
include("general_TMM.jl")
include("optics_functions.jl")
```

Add `Sheet` to the export list (after `Layer,` near line 9):

```julia
export Layer,
       Sheet,
       TransferResult,
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for the `Sheet constructors` testset (other testsets unchanged).

- [ ] **Step 6: Commit**

```bash
git add src/sheet.jl src/TransferMatrix.jl test/sheets.jl test/runtests.jl
git commit -m "feat: add Sheet type for 2D conductive sheets"
```

---

## Task 2: `sheet_matrix` — the interface matrix `G`

**Files:**
- Modify: `src/sheet.jl` (add `sheet_matrix`)
- Test: `test/sheets.jl`

- [ ] **Step 1: Write the failing test**

Append to `test/sheets.jl`:

```julia
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `UndefVarError: sheet_matrix not defined`.

- [ ] **Step 3: Add `sheet_matrix` to `src/sheet.jl`**

```julia
"""
    sheet_matrix(sheet::Sheet, λ)

Return the 4×4 interface matrix `G` for a conductive sheet at wavelength `λ`, in the
dynamical-matrix field basis `(Eₓ, Eᵧ, H_y, −Hₓ)`. This is the inverse of the physical
field-jump matrix, `G = G_phys⁻¹`, so its off-diagonal (H-row) entries are `+σ̃` with the
dimensionless conductivity `σ̃ = Z₀ σ` (`Z₀ = √(μ₀/ε₀)`). It is injected at an interface
as `Dᵢ⁻¹ G Dᵢ₊₁`. Tangential E (rows 1,2) is continuous; tangential H (rows 3,4) jumps.

`Z₀` is applied here and nowhere else (see [`Sheet`](@ref)).
"""
function sheet_matrix(sheet::Sheet, λ)
    σ = sheet.conductivity(λ)          # SMatrix{2,2,ComplexF64}, SI Siemens
    g = sqrt(μ_0 / ε_0) .* σ           # dimensionless; the ONLY Z₀ application
    z = zero(ComplexF64)
    o = one(ComplexF64)
    return @SMatrix [
        o       z       z  z;
        z       o       z  z;
        g[1,1]  g[1,2]  o  z;
        g[2,1]  g[2,2]  z  o
    ]
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for `sheet_matrix G structure`.

- [ ] **Step 5: Commit**

```bash
git add src/sheet.jl test/sheets.jl
git commit -m "feat: add sheet_matrix interface matrix G"
```

---

## Task 3: L–P reformulation + sheet injection in `transfer`

**Files:**
- Modify: `src/general_TMM.jl` (`_propagate_core`, `transfer`; add `_sheets_dict`, `_validate_sheet_indices`)
- Test: `test/sheets.jl`

This rewrites `_propagate_core` from the per-layer block form `Tᵢ = DᵢPᵢDᵢ⁻¹` to the
equivalent interface–propagation form so a sheet injects at one interface. The no-sheet
path stays numerically equivalent (guarded by the existing suite + the σ=0 test).

- [ ] **Step 1: Write the failing tests**

Append to `test/sheets.jl`:

```julia
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `transfer` has no `sheets` keyword (`MethodError`/unsupported kwarg).

- [ ] **Step 3: Add normalization + validation helpers to `src/general_TMM.jl`**

Add near the top of `src/general_TMM.jl` (after the struct definitions, before `poynting`):

```julia
# Normalize accepted sheet inputs (Dict or iterable of `i => sheet` pairs) to Dict{Int,Sheet}.
_sheets_dict(s::Dict{Int,Sheet}) = s
_sheets_dict(s) = Dict{Int,Sheet}(Int(k) => v for (k, v) in s)

# Unconditional structural validation: keys must index an interior interface.
function _validate_sheet_indices(sd, N)
    sd === nothing && return nothing
    for i in keys(sd)
        (1 ≤ i ≤ N - 1) || throw(ArgumentError("sheet index $i out of range; must be 1 ≤ i ≤ $(N - 1)"))
    end
    return nothing
end
```

- [ ] **Step 4: Rewrite `_propagate_core` to L–P form with sheet injection**

Replace the body of `_propagate_core` (currently `general_TMM.jl:289-318`) with:

```julia
function _propagate_core(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

    N = length(layers)
    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ε_0in = dielectric_constant(nx_in)
    ξ = √(ε_0in) * sin(θ)

    Λ_1324 = @SMatrix [1 0 0 0;
                       0 0 1 0;
                       0 1 0 0;
                       0 0 0 1]

    no_sheets = sheets === nothing || isempty(sheets)

    D_prev, _, γ_first, q_first = layer_matrices(layers[1], λ, ξ, μ)
    γ_last = γ_first
    q_last = q_first

    Γ = SMatrix{4,4,ComplexF64}(I)
    for i in 2:N
        layer = layers[i]
        D_cur, P_cur, γ_cur, q_cur = layer_matrices(layer, λ, ξ, μ)
        if no_sheets || !haskey(sheets, i - 1)
            L = D_prev \ D_cur                                  # interface (i-1, i)
        else
            L = D_prev \ (sheet_matrix(sheets[i - 1], λ) * D_cur)
        end
        Γ *= L                                                  # first ⇒ D₀⁻¹D₂ ; last ⇒ D_{N-1}⁻¹D_f
        if i < N
            Γ *= P_cur(layer.thickness)                         # propagate interior layer i
        end
        D_prev = D_cur
        if i == N
            γ_last = γ_cur
            q_last = q_cur
        end
    end

    Γ = (Λ_1324 \ Γ) * Λ_1324
    r, R, t, T = calculate_tr(Γ)
    S = poynting(ξ, q_first, q_last, γ_first, γ_last, t, r)

    return Γ, S
end
```

- [ ] **Step 5: Add `sheets` to `transfer`**

Modify `transfer`'s signature and the `_propagate_core` call (`general_TMM.jl:509-511`):

```julia
function transfer(λ, layers; θ=0.0, μ=1.0, sheets=nothing, validate::Bool=false)

    sd = sheets === nothing ? nothing : _sheets_dict(sheets)
    _validate_sheet_indices(sd, length(layers))
    Γ, S = _propagate_core(λ, layers; θ=θ, μ=μ, sheets=sd)
```

(The rest of `transfer` is unchanged for this task; `_validate_physics` gains sheet awareness in Task 8.)

- [ ] **Step 6: Run the tests to verify they pass**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for `sheet σ=0 regression` and `analytic conductive interface`, AND all
existing testsets still pass (the L–P refactor is behavior-preserving for no-sheet stacks).

**If the analytic test fails by a pure sign** (`Rss`/`Rpp` match the formula with `-g`
instead of `+g`): the injected `G` sign is inverted — flip the sign of the `g[...]`
entries in `sheet_matrix` (Task 2 Step 3) and re-run.

- [ ] **Step 7: Commit**

```bash
git add src/general_TMM.jl test/sheets.jl
git commit -m "feat: inject conductive sheets via L-P transfer reformulation"
```

---

## Task 4: Thread `sheets` through the sweeps

**Files:**
- Modify: `src/general_TMM.jl` (`_sweep_spectra`, `sweep_angle`, `sweep_thickness`)
- Test: `test/sheets.jl`

- [ ] **Step 1: Write the failing test**

Append to `test/sheets.jl`:

```julia
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `sweep_angle` has no `sheets` keyword.

- [ ] **Step 3: Add `sheets` to `_sweep_spectra` and forward it**

In `_sweep_spectra` (`general_TMM.jl:619`), add `sheets=nothing` to the signature and pass it
to the `transfer` call inside `compute_row` (`general_TMM.jl:638`):

```julia
function _sweep_spectra(outer_vals, inner_vals; threads::Bool=true, verbose::Bool=false, make_layers, angle_for, sheets=nothing)
```

```julia
            result = transfer(inner_vals[j], layers_i; θ=θ, sheets=sheets)
```

- [ ] **Step 4: Add `sheets` to `sweep_angle` and `sweep_thickness`**

`sweep_angle` (`general_TMM.jl:663`):

```julia
function sweep_angle(λs, θs, layers; sheets=nothing, threads::Bool=true, verbose::Bool=false)
    sd = sheets === nothing ? nothing : _sheets_dict(sheets)
    _validate_sheet_indices(sd, length(layers))
    return _sweep_spectra(θs, λs; threads=threads, verbose=verbose,
        make_layers = _ -> layers,
        angle_for = i -> θs[i],
        sheets = sd)
end
```

`sweep_thickness` (`general_TMM.jl:691`):

```julia
function sweep_thickness(λs, ts, layers, t_index::Int; θ=0.0, sheets=nothing, threads::Bool=true, verbose::Bool=false)
    sd = sheets === nothing ? nothing : _sheets_dict(sheets)
    _validate_sheet_indices(sd, length(layers))
    dispersion_func = layers[t_index].dispersion
    layers_base = collect(layers)

    return _sweep_spectra(ts, λs; threads=threads, verbose=verbose,
        make_layers = i -> begin
            layers_i = copy(layers_base)
            layers_i[t_index] = Layer(dispersion_func, ts[i])
            layers_i
        end,
        angle_for = _ -> θ,
        sheets = sd)
end
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for `sweep_angle forwards sheets`.

- [ ] **Step 6: Commit**

```bash
git add src/general_TMM.jl test/sheets.jl
git commit -m "feat: forward sheets through sweep_angle and sweep_thickness"
```

---

## Task 5: `_propagate_full` (with `qs`) and `propagate` wrapper

**Files:**
- Modify: `src/general_TMM.jl` (`propagate` → `_propagate_full` + wrapper)
- Test: `test/sheets.jl`

- [ ] **Step 1: Write the failing test**

Append to `test/sheets.jl`:

```julia
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `_propagate_full not defined`.

- [ ] **Step 3: Replace `propagate` with `_propagate_full` + wrapper**

Replace the entire `propagate` function (`general_TMM.jl:336-383`) with:

```julia
"""
    _propagate_full(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Internal full transfer-matrix pass. Returns `(Γ, S, Ds, Ps, γs, qs)` — like
[`propagate`](@ref) but also returns the per-layer eigenvalue vectors `qs`,
needed for magnetic-field reconstruction. Supports conductive sheets (Task: sheets).
"""
function _propagate_full(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

    N = length(layers)
    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ε_0in = dielectric_constant(nx_in)
    ξ = √(ε_0in) * sin(θ)

    Λ_1324 = @SMatrix [1 0 0 0;
                       0 0 1 0;
                       0 1 0 0;
                       0 0 0 1]

    no_sheets = sheets === nothing || isempty(sheets)

    D_1, P_1, γ_1, q_1 = layer_matrices(layers[1], λ, ξ, μ)
    Ds = Vector{typeof(D_1)}(undef, N)
    Ps = Vector{typeof(P_1)}(undef, N)
    γs = Vector{typeof(γ_1)}(undef, N)
    qs = Vector{typeof(q_1)}(undef, N)
    Ds[1] = D_1; Ps[1] = P_1; γs[1] = γ_1; qs[1] = q_1

    Γ = SMatrix{4,4,ComplexF64}(I)
    for i in 2:N
        D_i, P_i, γ_i, q_i = layer_matrices(layers[i], λ, ξ, μ)
        Ds[i] = D_i; Ps[i] = P_i; γs[i] = γ_i; qs[i] = q_i
        if no_sheets || !haskey(sheets, i - 1)
            L = Ds[i - 1] \ D_i
        else
            L = Ds[i - 1] \ (sheet_matrix(sheets[i - 1], λ) * D_i)
        end
        Γ *= L
        if i < N
            Γ *= P_i(layers[i].thickness)
        end
    end

    Γ = (Λ_1324 \ Γ) * Λ_1324
    r, R, t, T = calculate_tr(Γ)
    S = poynting(ξ, q_1, qs[N], γ_1, γs[N], t, r)

    return Γ, S, Ds, Ps, γs, qs
end

"""
    propagate(λ, layers; θ=0.0, μ=1.0, sheets=nothing)

Calculate the transfer matrix and Poynting vector for the structure, plus the
per-layer `D`, `P`, and `γ` matrices used for field reconstruction. Returns the
5-tuple `(Γ, S, Ds, Ps, γs)`. See [`transfer`](@ref) for the public R/T API.
"""
propagate(λ, layers; θ=0.0, μ=1.0, sheets=nothing) =
    _propagate_full(λ, layers; θ=θ, μ=μ, sheets=sheets)[1:5]
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for `_propagate_full and propagate contract`, and the existing
`integration.jl` tests that destructure `propagate` as a 5-tuple still pass.

- [ ] **Step 5: Commit**

```bash
git add src/general_TMM.jl test/sheets.jl
git commit -m "refactor: add _propagate_full with qs; keep propagate 5-tuple"
```

---

## Task 6: `_field` core; refactor `efield`; sheet support in fields

**Files:**
- Modify: `src/general_TMM.jl` (add `_field`; rewrite `efield` as a wrapper)
- Test: `test/sheets.jl`

- [ ] **Step 1: Write the failing tests**

Append to `test/sheets.jl`:

```julia
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
    dz = 0.001
    ef = efield(0.6, layers; dz = dz, sheets = sheets)

    # Find samples straddling the layer-2/3 interface (at z = boundaries[2])
    zint = ef.boundaries[2]
    jbelow = findlast(z -> z ≤ zint, ef.z)
    jabove = jbelow + 1
    # Tangential E (Ex, Ey) continuous across the sheet
    @test isapprox(ef.p[1, jbelow], ef.p[1, jabove]; atol = 5e-3)   # Ex (p)
    @test isapprox(ef.s[2, jbelow], ef.s[2, jabove]; atol = 5e-3)   # Ey (s)
end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `efield` has no `sheets` keyword yet.

- [ ] **Step 3: Add the `_field` core**

Add to `src/general_TMM.jl` (just before the current `efield` definition, ~line 739):

```julia
# Shared core for efield/hfield: runs _propagate_full once, performs the backward
# mode-coefficient recursion (with sheet injection), samples the z-grid, and returns
# everything both wrappers need. E and H differ only in the final per-z reconstruction.
function _field(λ, layers; θ=0.0, μ=1.0, dz=0.001, sheets=nothing)

    sd = sheets === nothing ? nothing : _sheets_dict(sheets)
    _validate_sheet_indices(sd, length(layers))
    no_sheets = sd === nothing || isempty(sd)

    Γ, S, Ds, Ps, γs, qs = _propagate_full(λ, layers; θ=θ, μ=μ, sheets=sd)
    r, R, t, T = calculate_tr(Γ)

    nx_in, _, _ = get_refractive_indices(layers[1], λ)
    ξ = √(dielectric_constant(nx_in)) * sin(θ)

    first_layer = layers[1]
    last_layer = layers[end]
    nlay = length(layers)

    Eplus_p = zeros(ComplexF64, nlay, 4)
    Eminus_p = zeros(ComplexF64, nlay, 4)
    Eplus_s = zeros(ComplexF64, nlay, 4)
    Eminus_s = zeros(ComplexF64, nlay, 4)

    Eplus_p[end, :] = [t[1], t[2], 0, 0]
    Eplus_s[end, :] = [t[3], t[4], 0, 0]

    P_f = Ps[end]
    Eminus_p[end, :] = P_f(last_layer.thickness) \ Eplus_p[end, :]
    Eminus_s[end, :] = P_f(last_layer.thickness) \ Eplus_s[end, :]

    D_i = Ds[end]
    for l in reverse(eachindex(layers))
        if l >= 2
            layer = layers[l - 1]
            P_prev = Ps[l - 1]
            if no_sheets || !haskey(sd, l - 1)
                L_i = Ds[l - 1] \ D_i
            else
                L_i = Ds[l - 1] \ (sheet_matrix(sd[l - 1], λ) * D_i)
            end
            Eminus_p[l - 1, :] = L_i * Eplus_p[l, :]
            Eminus_s[l - 1, :] = L_i * Eplus_s[l, :]
            Eplus_p[l - 1, :] = P_prev(layer.thickness) * Eminus_p[l - 1, :]
            Eplus_s[l - 1, :] = P_prev(layer.thickness) * Eminus_s[l - 1, :]
            D_i = Ds[l - 1]
        end
    end

    interface_positions, total_thickness = find_bounds(layers)
    interface_positions .-= first_layer.thickness
    zs = range(-first_layer.thickness, interface_positions[end], step=dz)

    nz = length(zs)
    amp_p = zeros(ComplexF64, 4, nz)
    amp_s = zeros(ComplexF64, 4, nz)
    layer_of_z = Vector{Int}(undef, nz)

    i = 1
    for (j, z) in enumerate(zs)
        if i < nlay && z > interface_positions[i]
            i += 1
        end
        P_i = Ps[i]
        amp_p[:, j] = P_i(-(z - interface_positions[i])) * Eminus_p[i, :]
        amp_s[:, j] = P_i(-(z - interface_positions[i])) * Eminus_s[i, :]
        layer_of_z[j] = i
    end

    return (; zs, boundaries = interface_positions[1:end - 1],
              amp_p, amp_s, layer_of_z, γs, qs, ξ, μ)
end
```

- [ ] **Step 4: Rewrite `efield` as a wrapper over `_field`**

Replace the body of `efield` (`general_TMM.jl:739-800`) with (keep its existing docstring,
add `sheets=nothing` to the signature):

```julia
function efield(λ, layers; θ=0.0, μ=1.0, dz=0.001, sheets=nothing)

    F = _field(λ, layers; θ=θ, μ=μ, dz=dz, sheets=sheets)
    nz = length(F.zs)
    p = zeros(ComplexF64, 3, nz)
    s = zeros(ComplexF64, 3, nz)

    for j in 1:nz
        γ = F.γs[F.layer_of_z[j]]
        ap = view(F.amp_p, :, j)
        as = view(F.amp_s, :, j)
        @views p[:, j] = ap[1] * γ[1, :] + ap[2] * γ[2, :] + ap[3] * γ[3, :] + ap[4] * γ[4, :]
        @views s[:, j] = as[1] * γ[1, :] + as[2] * γ[2, :] + as[3] * γ[3, :] + as[4] * γ[4, :]
    end

    return ElectricField(F.zs, p, s, F.boundaries)
end
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for both new testsets; the existing field tests in `integration.jl`
("single thin film propagate and field", "single thin film resonance and continuity")
still pass (efield output unchanged for no-sheet stacks).

- [ ] **Step 6: Commit**

```bash
git add src/general_TMM.jl test/sheets.jl
git commit -m "refactor: extract _field core; add sheet support to efield"
```

---

## Task 7: `MagneticField` and `hfield`

**Files:**
- Modify: `src/general_TMM.jl` (add `MagneticField`, `_h_eigvecs`, `hfield`)
- Modify: `src/TransferMatrix.jl` (export `hfield`, `MagneticField`)
- Test: `test/sheets.jl`

- [ ] **Step 1: Write the failing tests**

Append to `test/sheets.jl`:

```julia
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
    # R/T from transfer must equal R/T reconstructed from the field core's coefficients.
    air = Layer(λ -> complex(1.0), 0.0)
    spacer = Layer(λ -> complex(1.0), 0.5)
    sub = Layer(λ -> complex(1.5), 0.0)
    layers = [air, spacer, sub]
    sheets = Dict(2 => TransferMatrix.Sheet(2.0e-4 + 1.0e-4im))

    res = transfer(0.6, layers; θ = 0.1, sheets = sheets)
    Γ, _ = TransferMatrix._propagate_core(0.6, layers; θ = 0.1, sheets = TransferMatrix._sheets_dict(sheets))
    r, R, t, T = TransferMatrix.calculate_tr(Γ)
    @test isapprox(res.Rpp, R[1]; atol = 1e-12)
    @test isapprox(res.Rss, R[2]; atol = 1e-12)
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `hfield`/`MagneticField` not defined.

- [ ] **Step 3: Add `MagneticField`, `_h_eigvecs`, `hfield`**

Add `MagneticField` next to `ElectricField` in `src/general_TMM.jl` (after the
`ElectricField` struct, ~line 25):

```julia
"""
    MagneticField

Spatial magnetic-field profile through a layered structure, mirroring
[`ElectricField`](@ref). Fields are in impedance-normalized units `H̃ = Z₀ H_SI`
so `|E| ~ |H̃|` for a plane wave and E/H can be overlaid directly.

- `z`: position coordinates (same grid as `efield`)
- `p`: `(Hx, Hy, Hz)` for p-polarized incidence
- `s`: `(Hx, Hy, Hz)` for s-polarized incidence
- `boundaries`: z-positions of interfaces
"""
struct MagneticField{Z<:AbstractVector{Float64}}
    z::Z
    p::Matrix{ComplexF64}
    s::Matrix{ComplexF64}
    boundaries::Vector{Float64}

    function MagneticField(z::Z, p, s, boundaries) where {Z<:AbstractVector{Float64}}
        size(p, 2) == length(z) || throw(ArgumentError("p field columns must match z length"))
        size(s, 2) == length(z) || throw(ArgumentError("s field columns must match z length"))
        new{Z}(z, p, s, boundaries)
    end
end
```

Add the H-eigenvector helper and `hfield` after `efield` (~line 800):

```julia
# H eigenvectors per mode from the E eigenvectors γ and eigenvalues q:
# H_m = (1/μ)(-q γ₂, q γ₁ - ξ γ₃, ξ γ₂) = (Hx, Hy, Hz). Rows 2,1 match
# dynamical_matrix rows 3,4 (H_y and -Hx); row 3 (Hz) is (k×E)_z = ξ E_y.
function _h_eigvecs(γ, q, ξ, μ)
    η = @MMatrix zeros(ComplexF64, 4, 3)
    for m in 1:4
        η[m, 1] = (-q[m] * γ[m, 2]) / μ
        η[m, 2] = (q[m] * γ[m, 1] - ξ * γ[m, 3]) / μ
        η[m, 3] = (ξ * γ[m, 2]) / μ
    end
    return SMatrix(η)
end

"""
    hfield(λ, layers; θ=0.0, μ=1.0, dz=0.001, sheets=nothing)

Calculate the magnetic-field profile through the structure, returning a
[`MagneticField`](@ref). Shares the sampling grid with [`efield`](@ref) (same `dz`
and arguments), so E and H can be overlaid; the in-plane H discontinuity at a
conductive sheet equals the surface current `ẑ × (σ_s E∥)`.

# Units / normalization
H is returned in impedance-normalized units `H̃ = Z₀ H_SI` (`Z₀ = √(μ₀/ε₀)`), so
`|E| ~ |H̃|` for a plane wave. Arguments and conventions match [`efield`](@ref).
"""
function hfield(λ, layers; θ=0.0, μ=1.0, dz=0.001, sheets=nothing)

    F = _field(λ, layers; θ=θ, μ=μ, dz=dz, sheets=sheets)
    nz = length(F.zs)
    p = zeros(ComplexF64, 3, nz)
    s = zeros(ComplexF64, 3, nz)

    for j in 1:nz
        li = F.layer_of_z[j]
        η = _h_eigvecs(F.γs[li], F.qs[li], F.ξ, F.μ)
        ap = view(F.amp_p, :, j)
        as = view(F.amp_s, :, j)
        @views p[:, j] = ap[1] * η[1, :] + ap[2] * η[2, :] + ap[3] * η[3, :] + ap[4] * η[4, :]
        @views s[:, j] = as[1] * η[1, :] + as[2] * η[2, :] + as[3] * η[3, :] + as[4] * η[4, :]
    end

    return MagneticField(F.zs, p, s, F.boundaries)
end
```

- [ ] **Step 4: Export `hfield` and `MagneticField`**

In `src/TransferMatrix.jl`, add to the export list (near `ElectricField,` / `efield,`):

```julia
       ElectricField,
       MagneticField,
       efield,
       hfield,
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for `hfield basics and shared z-grid`, `transfer vs field cross-check`,
and `hfield H-jump across a sheet`.

- [ ] **Step 6: Commit**

```bash
git add src/general_TMM.jl src/TransferMatrix.jl test/sheets.jl
git commit -m "feat: add MagneticField and hfield with H-jump validation"
```

---

## Task 8: `validate=true` awareness of lossy sheets

**Files:**
- Modify: `src/general_TMM.jl` (`_validate_physics`, `transfer` call site)
- Test: `test/sheets.jl`

- [ ] **Step 1: Write the failing tests**

Append to `test/sheets.jl`:

```julia
@testset "validate=true with sheets" begin
    n0 = 1.0
    air = Layer(λ -> complex(n0), 0.0)
    spacer = Layer(λ -> complex(1.0), 0.3)
    sub = Layer(λ -> complex(n0), 0.0)
    layers = [air, spacer, sub]

    # Lossless sheet (purely imaginary σ): R + T = 1, no warning expected.
    lossless = Dict(2 => TransferMatrix.Sheet(0.0 + 2.0e-4im))
    res = transfer(0.6, layers; sheets = lossless, validate = true)
    @test isapprox(res.Rss + res.Tss, 1.0; atol = 1e-6)

    # Lossy sheet (Re σ > 0): R + T < 1 and must NOT raise/warn as an energy violation.
    lossy = Dict(2 => TransferMatrix.Sheet(2.0e-4 + 1.0e-4im))
    res2 = transfer(0.6, layers; sheets = lossy, validate = true)
    @test res2.Rss + res2.Tss < 1.0 + 1e-6
end
```

- [ ] **Step 2: Run the tests to verify they fail (or false-warn)**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: The lossy case currently triggers `_validate_physics` to treat the (lossless)
*layers* as lossless and warn "Energy conservation violated". The test asserting
`R + T < 1` passes, but to satisfy the intent (no false warning) we add sheet awareness.

- [ ] **Step 3: Make `_validate_physics` sheet-aware**

Change `_validate_physics`'s signature (`general_TMM.jl:544`) to accept `sheets`:

```julia
function _validate_physics(λ, layers, Tpp, Tss, Rpp, Rss; sheets=nothing, atol=1e-6, k_threshold=1e-10)
```

Replace the `is_lossless` computation (`general_TMM.jl:566-569`) with:

```julia
    layers_lossless = all(layers) do layer
        nx, ny, nz = get_refractive_indices(layer, λ)
        all(n -> abs(imag(n)) < k_threshold, (nx, ny, nz))
    end

    sheets_lossless = sheets === nothing || all(values(sheets)) do sheet
        σ = sheet.conductivity(λ)
        all(c -> abs(real(c)) < k_threshold, (σ[1,1], σ[1,2], σ[2,1], σ[2,2]))
    end

    is_lossless = layers_lossless && sheets_lossless
```

In `transfer`, forward the normalized `sd` to `_validate_physics` (`general_TMM.jl:525-527`):

```julia
    if validate
        _validate_physics(λ, layers, Tpp, Tss, Rpp, Rss; sheets=sd)
    end
```

- [ ] **Step 4: Run the tests to verify they pass (no false warning)**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for `validate=true with sheets`; no "Energy conservation violated"
warning for the lossy-sheet case.

- [ ] **Step 5: Commit**

```bash
git add src/general_TMM.jl test/sheets.jl
git commit -m "feat: make validate=true aware of lossy sheets"
```

---

## Task 9: Edge cases and thin-slab equivalence

**Files:**
- Test: `test/sheets.jl`

- [ ] **Step 1: Write the tests**

Append to `test/sheets.jl`:

```julia
@testset "sheet edge cases" begin
    n1, n2 = 1.0, 1.5
    a = Layer(λ -> complex(n1), 0.0)
    m1 = Layer(λ -> complex(1.4), 0.2)
    m2 = Layer(λ -> complex(1.6), 0.2)
    sub = Layer(λ -> complex(n2), 0.0)

    # Sheet at the first interface (i=1) and last (i=N-1) — touch semi-infinite layers
    layers = [a, m1, sub]
    r_first = transfer(0.6, layers; sheets = Dict(1 => TransferMatrix.Sheet(1.0e-4 + 0im)))
    r_last  = transfer(0.6, layers; sheets = Dict(2 => TransferMatrix.Sheet(1.0e-4 + 0im)))
    @test 0.0 ≤ r_first.Rss ≤ 1.0
    @test 0.0 ≤ r_last.Rss ≤ 1.0

    # Multiple sheets in one stack
    layers4 = [a, m1, m2, sub]
    r_multi = transfer(0.6, layers4; sheets = Dict(1 => TransferMatrix.Sheet(1.0e-4 + 0im),
                                                   3 => TransferMatrix.Sheet(2.0e-4 + 0im)))
    @test 0.0 ≤ r_multi.Rss ≤ 1.0
    @test 0.0 ≤ r_multi.Tss ≤ 1.0

    # Out-of-range index throws ArgumentError
    @test_throws ArgumentError transfer(0.6, layers; sheets = Dict(0 => TransferMatrix.Sheet(1.0e-4 + 0im)))
    @test_throws ArgumentError transfer(0.6, layers; sheets = Dict(3 => TransferMatrix.Sheet(1.0e-4 + 0im)))
end

@testset "anisotropic sheet cross-pol" begin
    a = Layer(λ -> complex(1.0), 0.0)
    sub = Layer(λ -> complex(1.5), 0.0)
    layers = [a, sub]

    diag_sheet = Dict(1 => TransferMatrix.Sheet(; xx = 2.0e-4, yy = 1.0e-4))
    res_d = transfer(1.0, layers; θ = π/6, sheets = diag_sheet)
    @test isapprox(res_d.Rps, 0.0; atol = 1e-10)
    @test isapprox(res_d.Rsp, 0.0; atol = 1e-10)

    offdiag_sheet = Dict(1 => TransferMatrix.Sheet(; xx = 2.0e-4, yy = 1.0e-4, xy = 1.0e-4, yx = 1.0e-4))
    res_o = transfer(1.0, layers; θ = π/6, sheets = offdiag_sheet)
    @test res_o.Rps > 1e-8 || res_o.Rsp > 1e-8
end

@testset "thin-slab equivalence" begin
    # Sheet(n,d) should match a thin Layer(n,d) as d -> 0.
    n = 2.0 + 0.1im
    a = Layer(λ -> complex(1.0), 0.0)
    sub = Layer(λ -> complex(1.0), 0.0)
    λ = 0.6
    θ = 0.1

    prev_err = Inf
    for d in (1e-3, 1e-4)
        slab = [a, Layer(λ0 -> n, d), sub]
        sheet = Dict(1 => TransferMatrix.Sheet(λ0 -> n, d))
        r_slab = transfer(λ, slab; θ = θ)
        r_sheet = transfer(λ, [a, sub]; θ = θ, sheets = sheet)
        err = abs(r_slab.Rss - r_sheet.Rss) + abs(r_slab.Tss - r_sheet.Tss)
        @test err < 10 * d            # error scales with d
        prev_err = err
    end
end
```

- [ ] **Step 2: Run the tests**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for `sheet edge cases`, `anisotropic sheet cross-pol`, `thin-slab equivalence`.
(If `thin-slab equivalence` is too tight, relax the bound to `< 50 * d` — the leading
correction is O(k₀d).)

- [ ] **Step 3: Commit**

```bash
git add test/sheets.jl
git commit -m "test: sheet edge cases, cross-pol, and thin-slab equivalence"
```

---

## Task 10: Docs autodocs wiring and full-suite check

**Files:**
- Modify: `docs/src/lib/public.md`, `docs/src/lib/internals.md`

- [ ] **Step 1: Add `sheet.jl` to the autodocs Pages**

In `docs/src/lib/public.md`, change the Transfer Matrix Functions `@autodocs` `Pages`:

```
Pages = ["general_TMM.jl", "layer.jl", "matrix_constructors.jl", "sheet.jl"]
```

In `docs/src/lib/internals.md`, change the `@autodocs` `Pages` the same way:

```
Pages = ["general_TMM.jl", "layer.jl", "matrix_constructors.jl", "sheet.jl"]
```

- [ ] **Step 2: Run the full test suite (includes Aqua)**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: ALL testsets PASS, including `Code quality (Aqua.jl)` (no method ambiguities,
no undefined exports, no stale deps from the new `Sheet`/`hfield`/`MagneticField` API).

If Aqua reports an ambiguity among `Sheet` constructors, narrow the offending method's
argument types (e.g. type `nx, ny` in `Sheet(nx, ny, d)` as `Union{Number,Function,RefractiveMaterial}`).

- [ ] **Step 3: Commit**

```bash
git add docs/src/lib/public.md docs/src/lib/internals.md
git commit -m "docs: register sheet.jl in API autodocs"
```

---

## Task 11: TMDC cavity example (non-CI)

**Files:**
- Create: `examples/tmdc_cavity.jl`

- [ ] **Step 1: Write the example**

Create `examples/tmdc_cavity.jl` (runs under the existing `examples/` environment with
CairoMakie; NOT added to the package test target, NOT a new root dependency):

```julia
# TMDC monolayer inside a Fabry–Pérot-like cavity, modeled as a conductive sheet.
# Demonstrates a Lorentzian-exciton sheet conductivity producing a polariton
# anticrossing in an angle sweep, plus an E/H field overlay.
#
# Run from the examples/ environment:
#   julia --project=examples examples/tmdc_cavity.jl

using TransferMatrix
using CairoMakie

# Lorentzian exciton sheet conductivity σ(λ) in SI Siemens.
# σ(ω) = i ε₀ c₀ f ω / (ω₀² - ω² - i ω Γ) with oscillator strength f (length units).
function exciton_sigma(λ)
    c0 = 299792458.0
    ε0 = 8.8541878128e-12
    λ0 = 0.62               # exciton wavelength (μm)
    Γλ = 0.01               # linewidth (μm)
    f = 5e-3                # oscillator strength (μm)
    ω  = 2π * c0 / λ
    ω0 = 2π * c0 / λ0
    Γ  = 2π * c0 * Γλ / λ0^2
    return im * ε0 * c0 * f * ω / (ω0^2 - ω^2 - im * ω * Γ)
end

# Simple half-wavelength cavity: two dielectric mirrors around a spacer, monolayer
# at the spacer center (a field antinode).
n_air, n_hi, n_spacer = 1.0, 2.2, 1.5
λ0 = 0.62
mirror_hi = Layer(λ -> complex(n_hi), λ0 / (4 * n_hi))
spacer_half = Layer(λ -> complex(n_spacer), λ0 / (2 * n_spacer) / 2)

air = Layer(λ -> complex(n_air), 0.0)
sub = Layer(λ -> complex(n_air), 0.0)
layers = [air, mirror_hi, spacer_half, spacer_half, mirror_hi, sub]
sheet_index = 3                          # interface between the two spacer halves (cavity center)
sheets = Dict(sheet_index => Sheet(exciton_sigma))

λs = range(0.58, 0.66, length = 400)
θs = range(0.0, 0.6, length = 200)
spec = sweep_angle(collect(λs), collect(θs), layers; sheets = sheets)

fig = Figure(size = (900, 400))
ax1 = Axis(fig[1, 1], xlabel = "wavelength (μm)", ylabel = "angle (rad)", title = "Rss (polariton anticrossing)")
heatmap!(ax1, λs, θs, spec.Rss', colormap = :viridis)

ef = efield(λ0, layers; dz = 0.001, sheets = sheets)
hf = hfield(λ0, layers; dz = 0.001, sheets = sheets)
ax2 = Axis(fig[1, 2], xlabel = "z (μm)", ylabel = "|field| (norm.)", title = "E / H at the monolayer")
lines!(ax2, ef.z, vec(sum(abs2, ef.s; dims = 1)) .^ 0.5, label = "|E| (s)")
lines!(ax2, hf.z, vec(sum(abs2, hf.s; dims = 1)) .^ 0.5, label = "|H| (s)")
for b in ef.boundaries
    vlines!(ax2, b, color = (:gray, 0.4))
end
axislegend(ax2)

save("tmdc_cavity.png", fig)
println("wrote tmdc_cavity.png")
```

- [ ] **Step 2: Smoke-check the example runs (manual, not CI)**

Run: `julia --project=examples examples/tmdc_cavity.jl`
Expected: writes `tmdc_cavity.png`; the heatmap shows an avoided crossing near λ₀.
(If the `examples/` env lacks deps, run `julia --project=examples -e 'using Pkg; Pkg.instantiate()'` first.)

- [ ] **Step 3: Commit**

```bash
git add examples/tmdc_cavity.jl
git commit -m "docs: add TMDC cavity (conductive sheet) example"
```

---

## Self-Review

**1. Spec coverage:**
- Sheet type + 4 constructor families + SI/Z₀ invariant → Task 1, Task 2 (Step 3 invariant). ✓
- `G = G_phys⁻¹` (+σ̃) interface matrix → Task 2. ✓
- L–P reformulation + sheet injection (transfer) → Task 3. ✓
- `sheets` normalization + unconditional index validation → Task 3 (helpers), used in Tasks 3,4,6. ✓
- Sweeps forward sheets (incl. `_sweep_spectra`) → Task 4. ✓
- `propagate` 5-tuple unchanged; `_propagate_full` adds qs → Task 5. ✓
- `_field` core; efield refactor + sheet support → Task 6. ✓
- `MagneticField` + `hfield` + H reconstruction + H-jump/cross-check/z-grid tests → Task 7. ✓
- `validate=true` lossy-sheet awareness → Task 8. ✓
- Analytic conductive-interface sign test, energy, thin-slab, cross-pol, edge cases → Tasks 3, 8, 9. ✓
- Docs autodocs wiring → Task 10. ✓
- TMDC cavity example (non-CI) → Task 11. ✓

**2. Placeholder scan:** Task 7 Step 1 deliberately shows a placeholder testset then immediately replaces it with the real one (so the file ends with the correct version) — acceptable as written; the executor pastes the replacement. No "TODO"/"implement later" remain.

**3. Type consistency:** `Sheet`, `sheet_matrix`, `_sheets_dict`, `_validate_sheet_indices`,
`_propagate_full`, `propagate` (5-tuple), `_field` (NamedTuple with `zs, boundaries, amp_p,
amp_s, layer_of_z, γs, qs, ξ, μ`), `_h_eigvecs`, `MagneticField`, `hfield`, `efield` —
names and signatures are consistent across tasks. `sheets` default sentinel is `nothing`
everywhere; normalization to `Dict{Int,Sheet}` happens once at each public boundary.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-06-05-conductive-sheets.md`. Two execution options:

1. **Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.
2. **Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints.

Which approach?
