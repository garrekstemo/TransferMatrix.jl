# TransferMatrix.jl

Berreman 4×4 transfer matrix method (TMM) for EM wave propagation in layered
optical media. R/T spectra, angle-resolved properties, E-field profiles,
anisotropic/birefringent media, 2D conductive sheets.

Entry point: `src/general_TMM.jl` (core algorithm). Browse `src/` for the rest
(`layer.jl`, `sheet.jl`, `matrix_constructors.jl`, `optics_functions.jl`).
General Julia/Makie/Pkg conventions live in the global `~/.claude/CLAUDE.md`.

## API

`transfer(λ, layers; θ=0.0, μ=1.0, sheets=nothing, validate=false)` returns a
`TransferResult{T}` struct with **8 fields**: `Tpp Tss Tps Tsp Rpp Rss Rps Rsp`
(includes cross-polarization terms). It does NOT return a tuple.

- `sweep_angle(λs, θs, layers; ...)` and
  `sweep_thickness(λs, ts, layers, t_index; θ=0.0, ...)` return a
  `TransferResult` whose fields are **matrices** (threaded by default).
- `efield(λ, layers; θ=0.0, μ=1.0, dz=0.001, sheets=nothing)` → `ElectricField`;
  `hfield` → `MagneticField`.
- Other exports: `Layer`, `Sheet`, `fresnel`, `fresnel_coefficients`, `airy`,
  `stopband`, `dbr_reflectivity`, `refractive_index`, `dielectric_constant`,
  `dielectric_tensor`, `find_bounds`, `isanisotropic`, `isrotated`,
  `euler_rotation_matrix`, `rotate_dielectric_tensor`, `get_refractive_indices`,
  `get_euler_angles`.

Examples and tutorials: `docs/src/guide/` (quickstart, tutorial, validation)
and `test/`.

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
conductivity tensor in Siemens**. Constructors: `Sheet(σ)`,
`Sheet(; xx, yy, xy=0, yx=0)`, `Sheet(material, d)`, `Sheet(nx, ny, d)`.
Pass via the `sheets=` kwarg (Dict or iterable of `i => sheet`) on
`transfer` / `sweep_angle` / `sweep_thickness` / `efield`. See `src/sheet.jl`.

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
  medium, so T needs a Poynting-vector ratio (`T = S_out,z / S_inc,z`).
  Reflectance uses `R = |r|²` directly.
- **γᵢ₃₃ sign correction**: Eq. 13 of Xu et al. (2000) / Eq. 20 of Passler (2017)
  has a sign error in γᵢ₃₃ as printed. All γⱼ₃ must satisfy the z-constraint
  `γⱼ₃ = −[(μεᵢ₃₁ + ξqⱼ)γⱼ₁ + μεᵢ₃₂γⱼ₂] / (μεᵢ₃₃ − ξ²)` (coefficient of γⱼ₂ is
  **minus**, not plus). Hard-won, not in the published papers — see `errata.md`.

Full equation reference (coordinate system, eigenmode sorting, γ/D/P matrices,
transfer construction, r/t coefficients, Fresnel, edge cases, energy
conservation): `.claude/rules/berreman-4x4-equations.md`.

## Known numerical limitations

- **#70** (resolved): rotated anisotropic crystals conserve energy. The quoted
  `Tpp+Rpp ≈ 0.998` (uniaxial `euler=(π/6,π/4,0)`) was a **budget artifact** — it
  omits the reflected cross-pol `Rps` that an out-of-plane tilt produces. The correct
  per-input-polarization budget is `Rpp+Rps+Tpp = 1` and `Rss+Rsp+Tss = 1` (holds to
  ~1e-14; the Poynting `Tpp` already includes the cross-*transmitted* power, so do
  **not** also add `Tps`/`Tsp`). Convention: `r_{in,out}` — p-input cross-reflection
  is `Rps`, not `Rsp` (they coincide only at normal incidence). Separately, a genuine
  bug — transmission into an **anisotropic substrate** (`Rpp+Rps+Tpp ≈ 1.017`) — was
  fixed in `poynting()`: the two transmitted substrate eigenmodes carry different
  wavevectors, so their Poynting vectors are now summed per-mode instead of using one
  wavevector for the combined field.
- **#71**: anisotropic ambient at oblique incidence → NaN (uses nx for ξ).
- **#72**: absorbing incident medium — |r|² is not a true energy reflectance. The
  Poynting vector is non-additive for absorbing incident media (interference
  cross-terms; Ortiz & Mochán 2005), so only the transmitted Poynting vectors are
  used for output.

## Issue-fixing workflow

Before fixing any issue, write out three steps: (1) the current situation/context,
(2) the problem and when/why it manifests, (3) the proposed fix and why it solves
the root cause. Log detailed fixes in `errata.md` with the same depth.

**`errata.md` is gitignored and local** — it is not checked in, so do not assume
it exists in a fresh checkout; create it if missing.
