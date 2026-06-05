# Design: 2D conductive sheets (TMDC monolayers) in the transfer matrix

**Date:** 2026-06-05
**Branch:** `feature/conductive-sheets` (off `origin/main`)
**Status:** Approved design → spec for implementation planning
**Revision:** Incorporates the 2026-06-05 multi-agent spec review (codebase-consistency,
Berreman/physics-correctness, completeness critic). Notably: the `G` sign is now
*derived* (not deferred); `propagate`'s public return is unchanged; `_sweep_spectra`
and `_validate_physics` wiring is made explicit; `_field`'s contract is specified.

## 1. Motivation

Atomically thin 2D materials (TMDC monolayers such as MoS₂/WS₂/WSe₂, graphene)
are ~0.6–0.7 nm thick. Modeling them as a finite dielectric slab forces an
arbitrary "effective thickness" paired with a bulk ε, where only the *product*
is physically meaningful at optical wavelengths. The natural description is an
infinitesimally thin sheet carrying a surface current **K** = σ_s **E**∥,
characterized by a **sheet (optical) conductivity** σ_s(ω). This is also how
2D-material optical constants are reported, and an excitonic resonance lives
naturally in σ_s(ω) — which produces the polariton anticrossing when a monolayer
is placed inside a cavity.

This feature adds conductive sheets to TransferMatrix.jl so a TMDC monolayer can
be simulated inside a Fabry–Pérot / DBR cavity, including reflection/transmission
spectra, angle/thickness sweeps, and electric **and** magnetic field profiles.

## 2. Physics background

### 2.1 Refractive index ↔ sheet conductivity

With the package's `exp(-iωt)` time convention, a slab of index `n` (ε = n²) and
thickness `d` is equivalent, in the `d → 0` limit, to a zero-thickness sheet with

```
σ_s(ω) = -i ω ε₀ d (n² - 1)          [SI, Siemens]      ω = 2π c₀ / λ
```

equivalently `ε(ω) = 1 + i σ_s / (ε₀ ω d)`. Lossy ⇒ `Re σ_s > 0` (consistent with
the package's convention that absorption gives `Im n > 0`, `Im ε > 0`).

### 2.2 The interface matrix `G` — sign derived, not deferred

A conducting sheet leaves tangential **E** continuous but makes tangential **H**
jump: `n̂ × (H⁺ − H⁻) = K = σ_s E∥`, with `n̂ = +ẑ` (incident → substrate). In the
package's **dynamical-matrix row basis `(Eₓ, Eᵧ, H_y, −Hₓ)`** — the order produced
by `dynamical_matrix` (rows 3,4 = `(qγ₁−ξγ₃)/μ = H_y` and `qγ₂/μ = −Hₓ`), *distinct*
from the Berreman eigenvector order `(Eₓ, H_y, Eᵧ, −Hₓ)` used in `poynting` — the
**physical field jump** `Ψ⁺ = G_phys Ψ⁻` is (with `σ̃ = Z₀ σ_s`, `Z₀ = √(μ₀/ε₀) ≈ 376.73 Ω`):

```
            ⎡  1       0      0   0 ⎤
G_phys  =   ⎢  0       1      0   0 ⎥        (−σ̃ in the H rows: H jumps by −K)
            ⎢ −σ̃ₓₓ   −σ̃ₓᵧ    1   0 ⎥
            ⎣ −σ̃ᵧₓ   −σ̃ᵧᵧ    0   1 ⎦
```

**The matrix injected into the transfer product is `G ≡ G_phys⁻¹`, not `G_phys`.**
Derivation: continuity of the physical field across the interface with a sheet is
`A_i c_i = G_phys A_{i-1} c_{i-1}` (the code's dynamical matrix obeys
`A_{i-1}c_{i-1} = A_i c_i` with **no** sheet). Solving, the factor that replaces the
interface block `D_{i-1}⁻¹ D_i` is `D_{i-1}⁻¹ G_phys⁻¹ D_i`. Since `G_phys` is
unit-triangular, `G_phys⁻¹` simply flips the off-diagonal sign:

```
        ⎡  1      0      0   0 ⎤
G  =    ⎢  0      1      0   0 ⎥     G = G_phys⁻¹  (the injected matrix; +σ̃)
        ⎢ +σ̃ₓₓ  +σ̃ₓᵧ    1   0 ⎥
        ⎣ +σ̃ᵧₓ  +σ̃ᵧᵧ    0   1 ⎦
```

The `Λ₁₃₂₄` reordering is an outer similarity transform applied identically with or
without sheets and cannot change this internal sign. **Conclusion: `s = +1`.** The
analytic test (§7 test 2) is a **regression guard that confirms** this sign, not the
thing that determines it. (Note: the energy test §7 test 4 is *sign-insensitive* —
both `±σ̃` conserve energy — so sign validation rests entirely on test 2.)

- Diagonal σ̃ₓₓ, σ̃ᵧᵧ act within the p-block (Eₓ↔H_y) and s-block (Eᵧ↔Hₓ): **no
  p–s mixing** for an isotropic sheet *between isotropic media*. (Anisotropic
  neighbor layers can still mix p/s via their `D`-matrices — see §7 test 7.)
- Off-diagonal σ̃ₓᵧ, σ̃ᵧₓ couple p and s: **polarization conversion** (populates the
  cross-pol outputs `Rps`, `Rsp`, `Tps`, `Tsp` the code already computes).
- `G` depends only on σ̃ — not on ξ or μ.

### 2.3 n,k → σ̃ in package units; the single-Z₀ invariant

`σ̃ = Z₀ σ_s = -i k₀ d (n² - 1)`, with `k₀ = ω/c₀ = 2π/λ`.

**Invariant:** `Sheet.conductivity(λ)` always returns **SI Siemens (σ_s)**; the
factor `Z₀` is applied **exactly once, in `sheet_matrix`** (§4). No other code path
multiplies by `Z₀`. This prevents a double-`Z₀` bug.

## 3. The `Sheet` type — `src/sheet.jl` (new file)

A new file keeps `Layer` focused (it is unchanged). Canonical internal form: a
callable `λ → SMatrix{2,2,ComplexF64}` returning the **SI sheet conductivity tensor
`[σxx σxy; σyx σyy]` in Siemens**.

```julia
struct Sheet{F}
    conductivity::F   # λ -> SMatrix{2,2,ComplexF64}, SI Siemens (see §2.3 invariant)
end
```

### 3.1 Constructors (precise disambiguation)

`Sheet(σ)` is **scalar-only** (isotropic). Tensors go through the keyword form.
This avoids fragile return-type probing.

```julia
Sheet(σ::Number)            # constant isotropic: wraps to λ -> SMatrix(σ,0,0,σ)
Sheet(σ::Function)          # σ(λ)->Number, isotropic: wraps to λ -> SMatrix(σ(λ),0,0,σ(λ))
Sheet(; xx, yy, xy=0, yx=0) # anisotropic tensor; each entry a Number OR λ->Number, each wrapped
Sheet(material::RefractiveMaterial, d::Real)  # n,k→σ (isotropic), d in μm
Sheet(n::Function, d::Real)                   # n(λ)->index, isotropic, d in μm
Sheet(nx, ny, d::Real)      # in-plane anisotropic n,k→σ: σxx from nx, σyy from ny; d in μm
```

- All paths construct the SI-Siemens `λ → SMatrix{2,2,ComplexF64}` form, mirroring
  the `refractive_index` factory pattern (`layer.jl:159-185`).
- Disambiguation is by arity + type: 1 positional (`Number`/`Function`) = scalar
  conductivity; ends in `d::Real` = index→σ; keywords = tensor. `Sheet(material, d)`
  dispatches on `material::RefractiveMaterial`; `Sheet(n::Function, d)` on the
  function path; both are 2-positional ending in `Real`, distinguished by the first
  arg's type.
- `d` is in μm (matching package length units).
- A units helper (σ in units of σ₀ = e²/4ℏ) and an `euler` kwarg for in-plane axis
  rotation are **deferred** to the interleaved-stack work (§9) — no gap in scope here.

### 3.2 Tests for constructors

Unit-test each constructor returns the expected `SMatrix{2,2,ComplexF64}` at a sample
λ; assert `Sheet(n, d)` agrees with the closed-form `σ̃ = -i k₀ d (n²-1)` (÷Z₀ back to
SI) at one λ (catches a sign/units or double-Z₀ error).

## 4. The `G` builder — `sheet_matrix` in `src/matrix_constructors.jl`

```julia
function sheet_matrix(sheet::Sheet, λ)
    σ  = sheet.conductivity(λ)              # 2×2, SI Siemens
    σ̃  = sqrt(μ_0 / ε_0) .* σ               # dimensionless; the ONLY Z₀ application
    # returns the 4×4 SMatrix G = G_phys⁻¹ (+σ̃ in rows 3,4) in basis (Eₓ,Eᵧ,H_y,−Hₓ)
end
```

Placed near `dynamical_matrix`. Returns `SMatrix{4,4,ComplexF64}`. No dependence on
ξ or μ. **Not exported** (consistent with `dynamical_matrix`/`propagation_matrix`).

## 5. TMM injection — `src/general_TMM.jl`

### 5.1 L–P reformulation (equivalent to current code)

Reformulate the core transfer product from per-layer blocks `Tᵢ = Dᵢ Pᵢ Dᵢ⁻¹` into
the mathematically **identical** interface–propagation (L–P) form, which makes
interfaces first-class:

```
Γ = L₁₂ P₂ L₂₃ P₃ … L_{N-1,N},     Lᵢ,ᵢ₊₁ = Dᵢ⁻¹ Dᵢ₊₁
```

**Equivalence:** the current code computes `D₀⁻¹ (∏_{i=2}^{N-1} Dᵢ Pᵢ Dᵢ⁻¹) D_f`
which expands to `(D₀⁻¹ D₂) P₂ (D₂⁻¹ D₃) P₃ … (D_{N-1}⁻¹ D_f)`. The streaming loop
below reproduces exactly this; the `D₀⁻¹` and `D_f` boundary factors are **absorbed**
into the first and last interface terms (no separate `D_0 \ (… * D_f)` step). A sheet
at interface (i, i+1) replaces `Lᵢ,ᵢ₊₁ → Dᵢ⁻¹ G Dᵢ₊₁` with `G = G_phys⁻¹` (§2.2).

`_propagate_core` (used by `transfer`, the hot path) keeps its **streaming,
no-per-layer-allocation** style; it retains `q` for the first and last layers for the
trailing `poynting` call:

```julia
# layers[1] = incident (D₀, semi-infinite); layers[N] = substrate (D_f, semi-infinite)
D_prev, _, γ_first, q_first = layer_matrices(layers[1], λ, ξ, μ)   # incident
local D_cur, γ_last, q_last
Γ = SMatrix{4,4,ComplexF64}(I)
for i in 2:N
    D_cur, P_cur, γ_cur, q_cur = layer_matrices(layers[i], λ, ξ, μ)
    L = D_prev \ D_cur                                # no-sheet interface  (i-1, i)
    if !isempty(sheets) && haskey(sheets, i-1)        # sheet at interface (i-1, i)
        L = D_prev \ (sheet_matrix(sheets[i-1], λ) * D_cur)
    end
    Γ *= L                                            # first iter ⇒ D₀⁻¹ D₂ ; last ⇒ D_{N-1}⁻¹ D_f
    i < N && (Γ *= P_cur(layers[i].thickness))        # propagate interior layer i (skip substrate)
    D_prev = D_cur
    i == N && (γ_last = γ_cur; q_last = q_cur)
end
Γ = (Λ_1324 \ Γ) * Λ_1324
# poynting(ξ, q_first, q_last, γ_first, γ_last, t, r) as today
```

**Requirements:**
- No-sheet path (`isempty(sheets)`) must be numerically identical to current output
  (regression test §7.1) and within perf parity (4×4 static solves).
- First and last layers remain semi-infinite (no `P`).

`propagate` (the allocating variant, used by `efield`/`_field` and by integration
tests) is refactored as follows to avoid breaking its callers:

- **`propagate` keeps its exact public 5-tuple return `(Γ, S, Ds, Ps, γs)`.**
  Its two test callers (`integration.jl:64`, `:175`) and any 5-tuple destructure are
  unaffected.
- A new internal `_propagate_full(λ, layers; θ, μ, sheets)` returns the 6-tuple
  `(Γ, S, Ds, Ps, γs, qs)` (adds per-layer `qs`, currently computed but discarded).
  `propagate` becomes a thin wrapper returning the first five. `_field` (§6) calls
  `_propagate_full` to obtain `qs`.

### 5.2 Placement keyword, normalization, validation

`sheets` maps **layer index `i` → `Sheet` at interface (i, i+1)** (the sheet on the
+z side of layer `i`). The public functions accept either `Dict(i => sheet)` or an
iterable of `Pair{Int,Sheet}` and **normalize to a canonical `Dict{Int,Sheet}` once
at the API boundary** (so internal code can rely on `haskey`). Default: empty
`Dict{Int,Sheet}()`; the `isempty` fast-path (§5.1) keeps the no-sheet hot path
allocation-free and numerically identical.

**Validation (unconditional, not tied to `validate`):** on entry, throw
`ArgumentError("sheet index $i out of range; must be 1 ≤ i ≤ $(N-1)")` for any key
outside `1 ≤ i ≤ N-1` (mirrors `Layer`'s unconditional thickness check,
`layer.jl:53`). The existing `validate` kwarg remains reserved for the physics
energy check.

New signatures (all backward compatible — `sheets` defaults to empty):

```julia
transfer(λ, layers; θ=0.0, μ=1.0, sheets=Dict{Int,Sheet}(), validate=false)
sweep_angle(λs, θs, layers; sheets=Dict{Int,Sheet}(), threads=true, verbose=false)
sweep_thickness(λs, ts, layers, idx; θ=0.0, sheets=Dict{Int,Sheet}(), threads=true, verbose=false)
efield / hfield (; …, sheets=Dict{Int,Sheet}())
propagate / _propagate_core / _propagate_full   # gain sheets kwarg
```

**`_sweep_spectra` wiring (was silently dropping sheets):** `sweep_angle` and
`sweep_thickness` both delegate to `_sweep_spectra`, which is the function that calls
`transfer` (`general_TMM.jl:638`). `_sweep_spectra` **must gain a `sheets` kwarg and
forward it** to the `transfer` call in `compute_row`. To avoid per-λ re-normalization
inside the threaded loop, normalize `sheets` to a `Dict` once before the loop and
pass the `Dict` down (re-normalizing a `Dict` is idempotent/cheap). `sheets` is
read-only and shared across threads — thread-safe like `layers`.

**`sweep_thickness` + sheets:** sheet keys are layer indices into the base `layers`
vector and are passed through untouched as the swept layer's thickness changes (its
index is stable). A sheet keyed at or adjacent to the swept layer index is physically
correct (its interface position moves) — covered by a test (§7).

### 5.3 `validate=true` with lossy sheets

`_validate_physics` currently determines `is_lossless` from `layers` only. A lossy
sheet on a lossless layer stack would (correctly) give R+T<1 but falsely trip the
"energy conservation violated" warning. Fix: pass `sheets` to `_validate_physics`;
treat the system as lossless only if **all layers AND all sheets** are lossless,
where a sheet is lossless ⟺ `abs(real(σ(λ))) < k_threshold` for every component.
Then R+T≈1 is asserted only when both hold; otherwise R+T≤1.

## 6. Field profiles — `efield` (E) and `hfield` (H)

### 6.1 Shared core `_field` — explicit contract

Factor a single internal `_field(λ, layers; θ, μ, dz, sheets)` that runs
`_propagate_full` (with `sheets`) once and performs the backward mode-coefficient
recursion + spatial sampling, returning everything both wrappers need:

```julia
# _field returns a NamedTuple:
(; zs,                 # Vector{Float64} sample positions (the SHARED grid)
   boundaries,         # interface positions (for the field structs)
   amp_p,              # 4 × length(zs) complex mode amplitudes at each z (p-incidence)
   amp_s,              # 4 × length(zs) complex mode amplitudes at each z (s-incidence)
   layer_of_z,         # length(zs) Int: which layer each z sample is in
   γs, qs,             # per-layer γ (4×3) and q (length-4) for reconstruction
   ξ, μ)               # for the H eigenvector
```

`efield` and `hfield` call `_field` and differ **only** in the final per-z
reconstruction (E uses γ; H uses the H eigenvector η below). Because both consume the
same `zs`, the grids are guaranteed identical (same array, not two `range` calls).

- `efield` output (`ElectricField{Z}` with `z, p, s, boundaries`) is **unchanged**:
  `E[:,j] = Σ_m amp_p[m,j] · γ_{layer_of_z[j]}[m,:]` (same expression as today).
- A regression test asserts current `efield` output is byte-identical after the
  `_field` extraction, independent of the σ=0 path (§7.1).

### 6.2 Sheet injection in the recursion (key convention)

The backward recursion forms the interface factor `Ds[l-1] \ D_i` between layer `l-1`
(−z side) and layer `l` (+z side) (`general_TMM.jl:767`). A sheet at that interface
has key `l-1` (matching the forward loop's `i-1` convention, §5.1):

```julia
G = (!isempty(sheets) && haskey(sheets, l-1)) ? sheet_matrix(sheets[l-1], λ) : I
L_i = Ds[l-1] \ (G * D_i)
```

**Invariant:** the resolved `G = G_phys⁻¹` is applied **identically** in
`_propagate_core`, `_propagate_full`, and `_field`. A cross-check test (§7) asserts
`transfer` R/T with a sheet equals R/T reconstructed from `_field` mode coefficients
with the same sheet — the only test that catches a transfer-vs-field placement/sign
mismatch.

### 6.3 H reconstruction (`hfield`, new)

Per mode `m`, the H eigenvector follows from `k × E` with `k = (ω/c)(ξ,0,q)`,
consistent with `dynamical_matrix` rows 3–4 (`H_y` and `−Hₓ`); the `H_z = ξγ₂/μ`
term is the `(k×E)_z = ξE_y` component:

```
H_m = (1/μ) · ( −q_m γ_{m,2},   q_m γ_{m,1} − ξ γ_{m,3},   ξ γ_{m,2} )
        =      (   Hₓ,              H_y,                       H_z      )
```

`H[:,j] = Σ_m amp_p[m,j] · η_{layer_of_z[j]}[m,:]` (and the s-incidence analogue),
where `η` is built from `γ, q, ξ, μ`. H is returned in impedance-normalized units
`H̃ = Z₀ H_SI`, so `|E| ~ |H̃|` for a plane wave and E/H overlay directly.

New symmetric API (non-breaking), mirroring `ElectricField` exactly:

```julia
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

hfield(λ, layers; θ=0.0, μ=1.0, dz=0.001, sheets=Dict{Int,Sheet}()) -> MagneticField
```

`hfield`'s docstring includes a "Units / normalization" section (matching `efield`)
stating `H̃ = Z₀ H_SI` and that E/H from `efield`/`hfield` share the z-grid and can be
overlaid. A test asserts `efield(args).z == hfield(args).z` for identical args.

### 6.4 Parallel workstream

The H-field track (`_propagate_full` qs, `_field` extraction, `MagneticField`,
`hfield`, H-jump + z-grid + transfer-vs-field tests) is independent of the R/T core
(§5) once the shared-state contract (qs ordering via `_propagate_full`; `sheets` key
convention) is fixed, and can be implemented/reviewed in parallel.

## 7. Tests — `test/sheets.jl` (registered in `test/runtests.jl`)

Pin concrete parameters/tolerances (no under-specified comparisons):

1. **σ = 0 regression + refactor regression:** (a) a `Sheet` with zero conductivity
   yields output identical (atol 1e-12) to the no-sheet path for `transfer`,
   `sweep_angle`, `efield`; (b) current `efield`/`transfer` output is unchanged after
   the `_field`/`_propagate_full` refactor *with no sheets at all*.
2. **Analytic conductive interface (confirms the `G` sign):** single sheet between
   media `n₁`/`n₂`, normal + oblique (e.g. θ=30°), s and p, vs closed form, atol 1e-8:
   - Normal: `t = 2n₁/(n₁+n₂+σ̃)`, `r = (n₁−n₂−σ̃)/(n₁+n₂+σ̃)`.
   - s-pol: `r_s = (n₁cosθ₁ − n₂cosθ₂ − σ̃)/(n₁cosθ₁ + n₂cosθ₂ + σ̃)`.
   - p-pol: `r_p = (n₂cosθ₁ − n₁cosθ₂ + σ̃cosθ₁cosθ₂)/(n₂cosθ₁ + n₁cosθ₂ + σ̃cosθ₁cosθ₂)`,
     using the package's `fresnel_coefficients` p-convention
     (`optics_functions.jl:93`: `rp=(n₂cosθᵢ−n₁cosθₜ)/(n₂cosθᵢ+n₁cosθₜ)`), into which
     `+σ̃cosθ₁cosθ₂` adds. (Verified energy-conserving by the physics review.)
3. **Thin-slab equivalence:** `Sheet(n, d)` vs a thin `Layer(n, d)`, `n=2+0.1i`,
   `d ∈ {1e-3, 1e-4}` μm, compare `Rpp/Tpp/Rss/Tss`; rtol tightens as `d→0`.
4. **Energy (sign-insensitive — not a sign test):** lossless sheet (purely imaginary
   σ) ⇒ assert R+T ≈ 1 (atol 1e-6 for isotropic-medium sheets); lossy sheet
   (Re σ > 0) ⇒ assert R+T < 1 strictly. Restrict to Re σ ≥ 0.
5. **Anisotropy / cross-pol:** off-diagonal σ_xy ⇒ nonzero `Rps`/`Rsp`; diagonal σ ⇒
   cross-pol ≈ 0.
6. **Field continuity / jump:** across the sheet, `efield` shows tangential **E**
   continuous; `hfield` shows the in-plane **H** jump `ΔH∥ = ẑ × (σ_s E∥)`
   (real-space validation). Also assert `efield(args).z == hfield(args).z`.
7. **Edge cases:** (a) sheet at first interface i=1 and last i=N-1 (touch the
   semi-infinite layers); (b) ≥2 sheets in one stack; (c) sheet adjacent to a
   zero-thickness layer (ambient layers built with thickness 0.0); (d) sheet at
   oblique incidence with an **anisotropic** neighbor (assert energy bound and
   expected cross-pol behavior); (e) out-of-range `i` ⇒ `ArgumentError`;
   (f) `validate=true` with a lossless sheet (R+T=1 path must not warn) and with a
   lossy sheet (must not falsely warn — §5.3).
8. **Transfer-vs-field cross-check:** R/T from `transfer` with a sheet equals R/T
   reconstructed from `_field` mode coefficients with the same sheet (catches
   placement/sign mismatch between the two code paths).

Non-CI **example** (`examples/tmdc_cavity.jl`): a Lorentzian exciton σ(λ) in a
DBR/Fabry–Pérot cavity → polariton anticrossing via `sweep_angle`, plus an E/H field
overlay showing the monolayer at an antinode. Runs under the existing `examples/`
environment (CairoMakie per the project Makie conventions), **excluded from
`Pkg.test()`**; **no new dependency added to the root `Project.toml`**.

## 8. Files & exports

| File | Change |
|---|---|
| `src/sheet.jl` | **new** — `Sheet` type + constructors + n,k→σ conversion |
| `src/matrix_constructors.jl` | add `sheet_matrix(sheet, λ)` (**not exported**) |
| `src/general_TMM.jl` | L–P reformulation; `_propagate_full` (qs) + `propagate` wrapper; `sheets` threading; `_sweep_spectra` `sheets` kwarg; `_validate_physics` sheet-loss awareness; `_field` core; `hfield`; `MagneticField` |
| `src/TransferMatrix.jl` | `include("sheet.jl")`; export `Sheet`, `hfield`, `MagneticField` (not `sheet_matrix`) |
| `docs/src/lib/public.md`, `internals.md` | add `"sheet.jl"` to `@autodocs` `Pages`; confirm `Sheet`/`hfield`/`MagneticField` docstrings render |
| `test/sheets.jl` | **new** — tests §7; register in `test/runtests.jl` |
| `examples/tmdc_cavity.jl` | **new** — demonstration (examples/ env, non-CI) |

## 9. Non-goals (this sketch)

- Interleaved-stack syntax (`StackElement` hierarchy) — deferred; the `Sheet` type and
  `G`-injection built here migrate cleanly to it later.
- `euler` kwarg for in-plane sheet axis rotation — deferred to the interleaved work.
- Units helper for σ in σ₀ = e²/4ℏ — optional follow-up.
- Absorbed-power / `power_entering` bookkeeping beyond existing R/T.
- Docs guide/tutorial section — optional follow-up (API autodocs are wired in §8).

## 10. Risks / open items

- **`G` sign:** resolved a priori (`s = +1`, injected `G = G_phys⁻¹`, §2.2); test §7.2
  is a regression guard. The energy test §7.4 cannot pin the sign (both signs conserve
  energy), so §7.2 is load-bearing for the sign.
- **L–P regression:** the no-sheet path must remain numerically identical; guarded by
  §7.1(a,b) plus existing integration tests.
- **`propagate` contract:** unchanged (5-tuple); `qs` flows only via `_propagate_full`.
  Verified against the two `integration.jl` callers.
- **`_sweep_spectra` / `_validate_physics` wiring:** explicit (§5.2, §5.3) — both were
  latent silent-failure points.
- **p-pol analytic convention:** pinned to the package's `fresnel_coefficients`
  convention (§7.2); physics review confirmed the form is correct and energy-conserving.
