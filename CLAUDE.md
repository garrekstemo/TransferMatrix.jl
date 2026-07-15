# TransferMatrix.jl

Berreman 4×4 transfer matrix method (TMM) for EM wave propagation in layered
optical media. R/T spectra, angle-resolved properties, E-field profiles,
anisotropic/birefringent media, 2D conductive sheets.

Entry point: `src/transfer.jl` (public R/T API); browse `src/` for the rest.

## API

`transfer(λ, layers; θ=0.0, μ=1.0, sheets=nothing, validate=false, method=:exp)` returns a
`TransferResult{T}` struct with **8 fields**: `Tpp Tss Tps Tsp Rpp Rss Rps Rsp`
(includes cross-polarization terms). It does NOT return a tuple.

- `sweep_angle(λs, θs, layers; ...)` and
  `sweep_thickness(λs, ts, layers, t_index; θ=0.0, ...)` return a
  `TransferResult` whose fields are **matrices** (threaded by default).
- `method=:exp` (default) propagates interior layers via the matrix exponential of
  the Berreman Δ matrix (degeneracy-immune, no eigenmode sorting); `method=:eig`
  selects the eigenmode path. Both agree to ~1e-12; boundaries always use the
  eigenmode path. See `src/matrix_constructors.jl` (`layer_transfer_exp`).
- Other exports (`Layer`, `Sheet`, `efield`/`hfield`, `fresnel`, …): see
  `src/TransferMatrix.jl` export list and docstrings.

### Per-layer permeability (`mu=`)

`Layer` accepts an optional `mu=`: `nothing` (default — use the global `μ=`
fallback passed to `transfer` / `sweep_*`), a scalar, a constant 3×3 tensor, or
a wavelength-dependent `λ -> 3×3 matrix`. A layer whose `mu` is not `nothing`
overrides the global value for that layer only. Tensor-μ layers (including a
magnetic substrate) are fully supported. A magnetic *ambient at oblique
incidence* shares the existing #71 k_par caveat: the conserved in-plane
wavevector k_par is computed from the ambient permittivity only.

Helpers `gyrotropic_tensor` and `polder_permeability`: see docstrings. Gotcha:
`polder_permeability` returns a function of **frequency** — wrap with a
λ-to-frequency conversion for the `mu=` kwarg (`TransferMatrix.c_0` is the
package's internal speed-of-light constant in **μm/s**; map wavelength →
frequency according to your own unit system).

### Renamed functions (deprecated → use)

| Old (deprecated) | New |
|------------------|-----|
| `calculate_tr`   | `transfer` |
| `electric_field` | `efield` |
| `angle_resolved` | `sweep_angle` |
| `tune_thickness` | `sweep_thickness` |

### 2D sheets (Sheet / `sheets=`)

`Sheet` is a zero-thickness 2D conductive layer (e.g. TMDC monolayer, graphene).
Internally a callable `λ -> SMatrix{2,2,ComplexF64}` returning the **SI sheet
conductivity tensor in Siemens**. Pass via the `sheets=` kwarg (Dict or iterable
of `i => sheet`) on `transfer` / `sweep_angle` / `sweep_thickness` / `efield`.
Constructors: see `src/sheet.jl`.

## Conventions and gotchas

- **Units**: micrometers for λ, layer thickness, and `dz` — they must match. `θ`
  in radians from the surface normal (RefractiveIndex.jl convention). Optional:
  `using Unitful` lets `λ`/thickness/`dz` carry units (and `λ` accepts
  wavenumber/frequency/energy, `θ` accepts angle units like `45u"°"`), normalized
  to μm/radians via the `UnitfulExt` extension.
- **Geometry**: light propagates +z; z=0 at the first interface. The **first and
  last layers are semi-infinite** — thickness ignored for propagation but must
  be > 0.
- **Time convention**: `exp(-iωt)` (Berreman / Passler & Paarmann), opposite of
  Yeh. Flips the sign in the propagation matrix — matters when comparing formulas.
- **Transmittance ≠ |t|²** in general: the transmitted wave is in a different
  medium, so T needs a Poynting-vector ratio (`T = S_out,z / S_inc,z`). This holds
  for **all four** channels: `Tpp/Tss/Tps/Tsp` are each the Poynting flux of one
  substrate eigenmode (evaluated with its own wavevector) over the incident flux.
  Reflectance uses `R = |r|²` directly.
- **γᵢ₃₃ sign correction**: Eq. 13 of Xu et al. (2000) / Eq. 20 of Passler (2017)
  has a sign error in γᵢ₃₃ as printed. All γⱼ₃ must satisfy the z-constraint
  `γⱼ₃ = −[(μεᵢ₃₁ + ξqⱼ)γⱼ₁ + μεᵢ₃₂γⱼ₂] / (μεᵢ₃₃ − ξ²)` (coefficient of γⱼ₂ is
  **minus**, not plus). Hard-won, not in the published papers — see `errata.md`.

Full equation reference (coordinate system, eigenmode sorting, γ/D/P matrices,
transfer construction, r/t coefficients, Fresnel, edge cases, energy
conservation): `.claude/rules/berreman-4x4-equations.md`.

## Energy accounting & known limitations

- **Energy budget is per input polarization**: `Rpp+Rps+Tpp+Tps = 1` and
  `Rss+Rsp+Tss+Tsp = 1` (~1e-14 for lossless stacks, including rotated
  anisotropic crystals). Every transmittance is a **per-mode Poynting flux** —
  each substrate eigenmode evaluated with its own wavevector. Don't quote
  `Tpp+Rpp` alone as an energy check; converting stacks put power in the cross
  channels. Convention: `r_{in,out}` — p-input cross-reflection is `Rps`, not
  `Rsp` (they coincide only at normal incidence).
- **#71**: anisotropic ambient at oblique incidence → NaN (uses nx for k_par).
- **#72**: absorbing incident medium — |r|² is not a true energy reflectance
  (Poynting non-additivity, interference cross-terms; Ortiz & Mochán 2005). Only
  transmitted Poynting vectors are used for output.
- **#107**: TIR into an **anisotropic substrate** with one propagating + one
  evanescent transmitted eigenmode → `calculate_q` throws "Mode sorting failed"
  (deliberately loud rather than a wrong mode count; pinned by regression test).
  Supporting evanescent transmitted modes in anisotropic media is the eventual fix.
- The default `method=:exp` propagator removed interior-layer mode-sorting
  failures (#92), but #71/#107 are **boundary** problems and remain; it also does
  not address cascade overflow (#88).

## Issue-fixing workflow

Before fixing any issue, write out three steps: (1) the current situation/context,
(2) the problem and when/why it manifests, (3) the proposed fix and why it solves
the root cause. Log detailed fixes in `errata.md` with the same depth.

**`errata.md` is gitignored and local** — it is not checked in, so do not assume
it exists in a fresh checkout; create it if missing.
