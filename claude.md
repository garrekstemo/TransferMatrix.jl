# TransferMatrix.jl

A Julia implementation of the 4×4 transfer matrix method (TMM) for electromagnetic wave propagation in layered optical media.

## What This Package Does

TransferMatrix.jl simulates how light (electromagnetic waves) propagates through stacks of thin films and layered materials. It calculates:

- **Reflection and transmission spectra** (R/T) as functions of wavelength
- **Angle-resolved optical properties** for arbitrary angles of incidence
- **Electric field spatial profiles** through the layered structure
- **Polarization-dependent behavior** for both s-polarized and p-polarized light

### Physical Phenomena Modeled

- Distributed Bragg Reflectors (DBRs) and photonic bandgaps
- Fabry-Pérot cavity resonances
- Thin film interference patterns
- Brewster angle effects
- Absorbing materials (complex refractive indices)
- Anisotropic/birefringent media

## Mathematical Approach

The implementation uses the **Berreman 4×4 matrix formalism**:

1. Construct dielectric tensors from material refractive indices
2. Build the Δ-matrix describing field propagation as a differential equation
3. Solve for eigenvalues (q-values) representing propagation modes
4. Compute dynamical matrices (D) at each interface
5. Construct propagation matrices (P) with phase factors across each layer
6. Multiply transfer matrices across all layers to get the total system matrix (Γ)
7. Extract reflection/transmission coefficients from Γ
8. Use Poynting vector calculations for accurate energy conservation

Key references: Berreman (1972), Passler & Paarmann (2017), Xu et al. (2000).

## Project Structure

```
src/
├── TransferMatrix.jl      # Main module, exports, physical constants
├── layer.jl               # Layer type, refractive index handling
├── optics_functions.jl    # Fresnel equations, DBR formulas, stopband calculations
├── matrix_constructors.jl # M, a, Δ, γ, D, P matrix construction functions
└── general_TMM.jl         # Core TMM algorithm: propagate, calculate_tr, sweeps, E-field

test/
├── runtests.jl            # Test runner
├── types.jl               # Layer type tests
├── layer.jl               # Layer utility function tests
├── functions.jl           # Internal function tests
├── optics_functions.jl    # Optics formula tests
└── integration.jl         # End-to-end TMM validation, energy conservation

docs/src/
├── index.md               # Overview and philosophy
├── guide/quickstart.md    # DBR example walkthrough
├── guide/tutorial.md      # Brewster angle, Fabry-Pérot, E-field tutorials
├── lib/public.md          # API reference
└── bibliography.md        # Academic references
```

## Key Types

### `Layer{F,T}`
Represents a single optical layer with a dispersion function `F` (refractive index vs wavelength) and thickness `T`. Create with `Layer(material, thickness)` where material can be a `RefractiveMaterial` from RefractiveIndex.jl or a custom function.

### `Spectra`
Holds computed reflection/transmission values: `Rpp`, `Rss`, `Tpp`, `Tss` (p and s polarizations).

### `ElectricField`
Contains spatial E-field distribution: z-positions, field components (Ex, Ey, Ez) for each polarization, and layer boundary positions.

### `PropagationMatrix{Tω,Tq}`
Internal callable type that applies phase propagation over distance z.

## Main API Functions

| Function | Purpose |
|----------|---------|
| `Layer(material, thickness)` | Create a layer from material dispersion and thickness |
| `calculate_tr(λ, layers; θ=0.0, μ=1.0)` | Compute R/T at specific wavelength and angle |
| `sweep_angle(λs, θs, layers)` | Compute R/T over wavelength and angle arrays |
| `sweep_thickness(λs, ts, layers, idx; θ=0.0)` | Sweep layer thickness and compute R/T |
| `electric_field(λ, layers; θ=0.0, μ=1.0, dz=0.001)` | Compute E-field spatial profile |
| `refractive_index(material)` | Convert RefractiveMaterial to dispersion function |
| `dielectric_constant(n)` | Convert refractive index to dielectric constant |
| `dielectric_tensor(ε1, ε2, ε3)` | Create diagonal dielectric tensor |
| `find_bounds(layers)` | Get z-positions of layer interfaces |
| `fresnel(θ, n1, n2)` | Fresnel reflectance at an interface |
| `stopband(n1, n2)` | DBR photonic stopband bandwidth |
| `dbr_reflectivity(no, ns, n1, n2, N)` | Approximate DBR reflectivity |

## Usage Example

```julia
using TransferMatrix
using RefractiveIndex

# Define materials
n_air = RefractiveMaterial("other", "air", "Ciddor")
n_tio2 = RefractiveMaterial("main", "TiO2", "Sarkar")
n_sio2 = RefractiveMaterial("main", "SiO2", "Malitson")

# Build layer stack (DBR example)
λ_center = 0.633  # μm
air = Layer(n_air, 0.1)
tio2 = Layer(n_tio2, λ_center / (4 * 2.23))
sio2 = Layer(n_sio2, λ_center / (4 * 1.46))
substrate = Layer(n_sio2, 0.5)

layers = [air, tio2, sio2, tio2, sio2, tio2, sio2, substrate]

# Calculate at single wavelength
Tpp, Tss, Rpp, Rss = calculate_tr(0.633, layers, θ=0.0)

# Spectral sweep
λs = range(0.4, 0.9, length=500)
θs = [0.0]
spectra = sweep_angle(λs, θs, layers)

# Electric field profile
field = electric_field(0.633, layers, θ=0.0, dz=0.001)
```

## Units and Conventions

### Units
All length quantities must use **consistent units**. Following RefractiveIndex.jl conventions, **micrometers (μm)** are recommended:

| Quantity | Units | Notes |
|----------|-------|-------|
| Wavelength (`λ`) | μm | Must match layer thickness units |
| Layer thickness | μm | Must match wavelength units |
| Spatial step (`dz`) | μm | For electric field sampling |
| Angle (`θ`) | radians | Measured from surface normal |
| Transmittance/Reflectance | dimensionless | Values between 0 and 1 |

### Wave Propagation Convention

```
          incident medium (layer 1)
               ↓ light
    ─────────────────────────── z = 0 (first interface)
          layer 2
    ───────────────────────────
          layer 3
    ───────────────────────────
               ⋮
    ───────────────────────────
          substrate (last layer)
               ↓ transmitted light
```

- **Light propagates in the +z direction** (from first layer toward last layer)
- The **first and last layers are semi-infinite** (their thickness values are ignored for propagation, but thickness > 0 is required)
- **z = 0** is at the first interface (between layer 1 and layer 2)
- **Negative z** values are inside the incident medium
- **θ** is measured from the surface normal (z-axis); θ = 0 is normal incidence

### Time Convention and Propagation Matrix

The code uses the **exp(-iωt)** time convention, consistent with Berreman (1972) and Passler & Paarmann (2017).

The propagation matrix uses:
```
P(z) = diag(exp(-iωq₁z/c), exp(-iωq₂z/c), exp(-iωq₃z/c), exp(-iωq₄z/c))
```

This means:
- **Transmitted modes** (q with positive real part or positive imaginary part): field decays/propagates in +z direction
- **Reflected modes** (q with negative real part or negative imaginary part): field decays/propagates in -z direction

Note: Some references (e.g., Yeh) use exp(+iωt) time convention, which flips the sign in the propagation matrix. The physics is equivalent but care must be taken when comparing formulas.

### Polarization Convention
- **p-polarized (TM)**: Electric field in the plane of incidence (x-z plane)
- **s-polarized (TE)**: Electric field perpendicular to the plane of incidence (y direction)

## Physical Constants

Defined in `TransferMatrix.jl`:
- `ε_0 = 8.8541878128e-12` F/m (permittivity of free space)
- `μ_0 = 1.25663706212e-6` H/m (permeability of free space)
- `c_0 = 299792458` m/s (speed of light)

## Dependencies

- **RefractiveIndex.jl**: Material database for refractive indices
- **DataInterpolations.jl**: Interpolation for custom dispersion data
- **StaticArrays.jl**: High-performance fixed-size arrays
- **LinearAlgebra**: Matrix operations and eigendecomposition

## Internal Architecture

The TMM calculation pipeline in `general_TMM.jl`:

1. `calculate_tr` orchestrates the full calculation
2. `propagate` multiplies transfer matrices across the layer stack
3. Matrix constructors build intermediate matrices:
   - `Mmatrix`: 6×6 from dielectric/permeability tensors
   - `amatrix`: 6×6 intermediate matrix
   - `Δmatrix`: 4×4 reordered propagation matrix
   - `γmatrix`: 4×3 electric field components (with singularity handling)
   - `Dmatrix`: 4×4 dynamical matrix at interfaces
   - `Pmatrix`: 4×4 diagonal phase propagation
4. `poynting` calculates energy flow for accurate transmittance

## Testing

Run tests with:
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Key test areas:
- Layer construction and interface finding
- Fresnel equation validation against analytic formulas
- Energy conservation (R + T ≈ 1 for non-absorbing media)
- Field continuity at interfaces
- DBR stopband calculations

## Development Notes

### Issue-Fixing Workflow

Before fixing any issue in the codebase, follow this three-step explanation process:

1. **Explain the current situation**: Describe what the code currently does and the context around it.
2. **Explain the problem and why it's a problem**: Identify what's wrong, when it manifests, and what consequences it has (crashes, incorrect results, performance issues, etc.).
3. **Explain the fix and why it solves the problem**: Describe the proposed solution and how it addresses the root cause.

This ensures understanding before modification and creates documentation for future reference. Document fixes in ISSUES.md with the same level of detail.

### Technical Notes

- The code handles numerical singularities (degenerate q-values) using the Xu et al. (2000) approach
- Birefringence is detected automatically by comparing Poynting vector ratios
- Mode sorting uses thresholds on real/imaginary parts to distinguish transmitted vs reflected modes
- `sweep_angle` and `sweep_thickness` are designed to be thread-safe for parallel computation

## Mathematical Formalism (Passler & Paarmann 2017)

This section documents the key equations from the reference papers for debugging and development.

### Coordinate System & Conventions

- **Coordinates**: Multilayer surfaces parallel to x-y plane, z points from incident medium → substrate
- **Layer indexing**: i = 0 (incident), i = 1...N (layers), i = N+1 (substrate)
- **Wave vector**: kᵢ = (ω/c)(ξ, 0, qᵢ) where ξ = √ε_inc sin(θ) is conserved across all layers
- **Field vector ordering**: Ψ = (Eₓ, Hᵧ, Eᵧ, -Hₓ)ᵀ

### Eigenmode Sorting (Critical for Numerical Stability)

The four eigenvalues qᵢⱼ (j=1,2,3,4) of the Δ-matrix must be sorted consistently:

**Transmitted vs Reflected (Eq. 12)**:
- If qᵢⱼ is real: qᵢⱼ ≥ 0 → transmitted (j=1,2), qᵢⱼ < 0 → reflected (j=3,4)
- If qᵢⱼ is complex: Im(qᵢⱼ) ≥ 0 → transmitted, Im(qᵢⱼ) < 0 → reflected

**p vs s polarization sorting** using ratio C(qᵢⱼ):
```
C(qᵢⱼ) = |Ψᵢⱼ₁|² / (|Ψᵢⱼ₁|² + |Ψᵢⱼ₃|²)
```
Sort such that: C(qᵢ₁) > C(qᵢ₂) and C(qᵢ₃) > C(qᵢ₄)

Result: q₁ = p-trans, q₂ = s-trans, q₃ = p-refl, q₄ = s-refl

**For birefringent media**, use Poynting vector components instead:
```
C(qᵢⱼ) = |Sᵢⱼₓ|² / (|Sᵢⱼₓ|² + |Sᵢⱼᵧ|²)
```

### γ-matrix (Electric Field Eigenvectors)

For singularity-free handling (Xu et al. 2000), the γ components have conditional forms:

**Fixed components**: γᵢ₁₁ = γᵢ₂₂ = γᵢ₄₂ = -γᵢ₃₁ = 1

**Degenerate case** (qᵢ₁ = qᵢ₂ for isotropic/aligned materials):
```
γᵢ₁₂ = 0
γᵢ₁₃ = -(μᵢεᵢ₃₁ + ξqᵢ₁)/(μᵢεᵢ₃₃ - ξ²)
```

**Non-degenerate case** (qᵢ₁ ≠ qᵢ₂):
```
γᵢ₁₂ = [μᵢεᵢ₂₃(μᵢεᵢ₃₁ + ξqᵢ₁) - μᵢεᵢ₂₁(μᵢεᵢ₃₃ - ξ²)] / [(μᵢεᵢ₃₃ - ξ²)(μᵢεᵢ₂₂ - ξ² - qᵢ₁²) - μᵢ²εᵢ₂₃εᵢ₃₂]
γᵢ₁₃ = -(μᵢεᵢ₃₁ + ξqᵢ₁)/(μᵢεᵢ₃₃ - ξ²) - [μᵢεᵢ₃₂/(μᵢεᵢ₃₃ - ξ²)]γᵢ₁₂
```

Similar formulas for γᵢ₂₁, γᵢ₂₃, γᵢ₃₂, γᵢ₃₃, γᵢ₄₁, γᵢ₄₃ (see Eq. 20 in paper).

**Normalization (from 2019 erratum)**: γ_hat_ij = γ_vec_ij / |γ_vec_ij|

### D-matrix (Dynamical Matrix, 4×4)

```
       ⎛ γᵢ₁₁    γᵢ₂₁    γᵢ₃₁    γᵢ₄₁  ⎞
Aᵢ =   ⎜ γᵢ₁₂    γᵢ₂₂    γᵢ₃₂    γᵢ₄₂  ⎟
       ⎜ (qᵢ₁γᵢ₁₁-ξγᵢ₁₃)/μᵢ  ...        ⎟
       ⎝ qᵢ₁γᵢ₁₂/μᵢ    ...              ⎠
```

Row 3: (qᵢⱼγᵢⱼ₁ - ξγᵢⱼ₃)/μᵢ for each column j
Row 4: qᵢⱼγᵢⱼ₂/μᵢ for each column j

### P-matrix (Propagation Matrix, 4×4 diagonal)

```
Pᵢ = diag(exp(-iωqᵢ₁dᵢ/c), exp(-iωqᵢ₂dᵢ/c), exp(-iωqᵢ₃dᵢ/c), exp(-iωqᵢ₄dᵢ/c))
```

### Transfer Matrix Construction

- **Interface matrix**: Lᵢ = Aᵢ₋₁⁻¹ Aᵢ
- **Single layer**: Tᵢ = Aᵢ Pᵢ Aᵢ⁻¹
- **Full system**: Γ_N = A₀⁻¹ T_tot A_{N+1} = L₁P₁L₂P₂...L_NP_NL_{N+1}

**Reordering for Yeh convention**: Γ*_N = Λ₁₃₂₄⁻¹ Γ_N Λ₁₃₂₄

where Λ₁₃₂₄ swaps field ordering from (p-trans, s-trans, p-refl, s-refl) to (p-trans, p-refl, s-trans, s-refl)

### Reflection/Transmission Coefficients

From Γ* matrix elements (denominator D = Γ*₁₁Γ*₃₃ - Γ*₁₃Γ*₃₁):

```
rₚₚ = (Γ*₂₁Γ*₃₃ - Γ*₂₃Γ*₃₁) / D     tₚₚ = Γ*₃₃ / D
rₛₛ = (Γ*₁₁Γ*₄₃ - Γ*₄₁Γ*₁₃) / D     tₛₛ = Γ*₁₁ / D
rₚₛ = (Γ*₄₁Γ*₃₃ - Γ*₄₃Γ*₃₁) / D     tₚₛ = -Γ*₃₁ / D
rₛₚ = (Γ*₁₁Γ*₂₃ - Γ*₂₁Γ*₁₃) / D     tₛₚ = -Γ*₁₃ / D
```

**Reflectivity**: Rₖₗ = |rₖₗ|²

**Important**: Tₖₗ ≠ |tₖₗ|² in general! True transmittance requires Poynting vector calculation.

### Fresnel Equations (Single Interface)

For interface between media with refractive indices n₁ and n₂ at angle θ₁:

```
cos(θ₂) = √(1 - (n₁/n₂)²sin²(θ₁))   [Snell's law]

rₛ = (n₁cos(θ₁) - n₂cos(θ₂)) / (n₁cos(θ₁) + n₂cos(θ₂))
rₚ = (n₂cos(θ₁) - n₁cos(θ₂)) / (n₂cos(θ₁) + n₁cos(θ₂))

Rₛ = |rₛ|²
Rₚ = |rₚ|²
```

### Physical Edge Cases

1. **Grazing incidence** (θ → 90°): ξ → n₀, transmitted modes become evanescent
2. **Total internal reflection**: When sin(θ) > n₂/n₁, cos(θ₂) becomes imaginary, q purely imaginary
3. **Normal incidence** (θ = 0): ξ = 0, p and s polarizations degenerate, Rₚ = Rₛ
4. **Brewster angle**: rₚ = 0 when n₁cos(θ₂) = n₂cos(θ₁), i.e., tan(θ_B) = n₂/n₁
5. **Degenerate q-values**: Isotropic materials have qᵢ₁ = qᵢ₂ and qᵢ₃ = qᵢ₄ (requires special γ formulas)
6. **Critical angle**: θ_c = arcsin(n₂/n₁) for n₁ > n₂

### Energy Conservation

For non-absorbing media: R + T = 1

The true transmittance requires accounting for the change in wave impedance:
```
T = (n₂cos(θ₂))/(n₁cos(θ₁)) × |t|²
```

This is handled via Poynting vector calculations in the code.
