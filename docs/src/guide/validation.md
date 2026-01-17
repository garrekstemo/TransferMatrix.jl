# Physics Validation

TransferMatrix.jl includes built-in physics validation to help catch numerical issues and verify that results are physically meaningful. This page explains what gets validated and how to interpret warnings.

## Enabling Validation

Enable validation by passing `validate=true` to `transfer`:

```julia
result = transfer(λ, layers; validate=true)
```

Validation is disabled by default to avoid performance overhead in production calculations.

## Physical Constraints

The validation checks four physical constraints:

### NaN Detection

Detects whether any of `Tpp`, `Tss`, `Rpp`, `Rss` are `NaN`. `NaN` values indicate a numerical failure somewhere in the calculation pipeline. This can occur when:
- Layer parameters create a singular or near-singular transfer matrix
- Extreme refractive index values or layer thicknesses cause numerical overflow

**Example warning:**
```
┌ Warning: NaN detected in R/T values
│   Tpp = NaN
│   Tss = 0.85
│   Rpp = NaN
│   Rss = 0.15
└ @ TransferMatrix
```

### Bounds Check: 0 ≤ R, T ≤ 1

Reflectance and transmittance must be between 0 and 1. These are ratios of energy flux. Negative values indicate a sign error in the Poynting vector calculation. Values greater than 1 indicate numerical instability or an error in the transfer matrix construction.

**Common causes of violations:**
- Sign convention errors in the formalism
- Incorrect mode sorting (transmitted vs. reflected modes swapped)
- Numerical overflow in the propagation matrix for very thick layers

**Example warning:**
```
┌ Warning: Transmittance Tpp out of bounds [0, 1]
│   Tpp = -0.02
└ @ TransferMatrix
```

### Energy Conservation: R + T = 1 (Lossless Media)

For non-absorbing media, the sum of reflectance and transmittance must equal 1.
This is a fundamental consequence of energy conservation. In a lossless system, all incident power must either be reflected or transmitted—none is absorbed.
**The physics:** The Poynting vector ``\mathbf{S} = \frac{1}{2}\text{Re}(\mathbf{E} \times \mathbf{H}^*)`` represents the time-averaged power flow. The z-component gives the power flux through a surface. For energy conservation:

```math
S_{\text{incident}} = S_{\text{reflected}} + S_{\text{transmitted}}
```

Dividing by ``S_{\text{incident}}``:

```math
1 = R + T
```

**Important note:** The *amplitude* transmission coefficient ``|t|^2`` is **not** equal to the transmittance ``T`` in general. The true transmittance must account for the change in wave impedance between media:

```math
T = \frac{n_2 \cos\theta_2}{n_1 \cos\theta_1} |t|^2
```

This is why TransferMatrix.jl uses Poynting vector calculations rather than simply squaring the transmission coefficient.

**Example warning:**
```
┌ Warning: Energy conservation violated for p-polarization
│   Tpp = 0.72
│   Rpp = 0.25
│   sum = 0.97
│   expected = 1.0
│   deviation = 0.03
└ @ TransferMatrix
```

### Absorption Bound: R + T ≤ 1 (Lossy Media)

For absorbing media (where any layer has a non-zero extinction coefficient), the sum R + T must not exceed 1.
In absorbing media, some energy is converted to heat. The absorbed fraction is:

```math
A = 1 - R - T
```

If R + T > 1, it would imply negative absorption (energy creation), which is unphysical.

**Note:** The validation does *not* check that R + T < 1 for lossy media—only that it doesn't exceed 1. A value of exactly 1.0 could occur if the absorbing layer is very thin or if the imaginary part of the refractive index is negligible.

**Example warning:**
```
┌ Warning: Absorption violation: R + T > 1 for p-polarization
│   Tpp = 0.65
│   Rpp = 0.42
│   sum = 1.07
└ @ TransferMatrix
```

## When to Use Validation

**Use `validate=true` when:**
- Developing new simulations or testing new layer configurations
- Debugging unexpected results
- Working with edge cases (grazing incidence, very thin layers, high refractive index contrast)
- Validating against experimental data or other simulation tools

**Skip validation (`validate=false`, the default) when:**
- Running production calculations where you've already verified the setup
- Performance is critical (e.g., large parameter sweeps)
- You're certain the layer stack is well-behaved

## Interpreting Validation Failures

If you see validation warnings, here are some debugging steps:

1. **Check layer definitions:** Ensure thicknesses are positive and refractive indices are physically reasonable.

2. **Check units:** All lengths (wavelength, thickness) must use the same units. Micrometers (μm) are recommended.

3. **Check for edge cases:**
   - Grazing incidence (θ → 90°) can cause numerical issues
   - Very thin layers (thickness ≪ wavelength) may require higher precision
   - Very thick layers (thickness ≫ wavelength) can cause overflow in the propagation matrix

4. **Simplify the problem:** Test with a single interface (two layers) to verify basic functionality.

5. **Compare methods:** The code calculates R from transfer matrix coefficients and T from Poynting vectors. Significant disagreement suggests an issue with one of the calculation paths.

## Limitations

Validation does **not** catch all possible errors:

- **Field continuity:** Tangential E and H fields should be continuous across interfaces. This is tested in the test suite but not validated at runtime.
- **Incorrect but consistent results:** If there's a systematic error that affects both R and T equally, energy conservation might still hold.
- **Material data errors:** Validation can't detect if your refractive index data is wrong—only that the calculation is self-consistent.
- **Anisotropic media issues:** The validation assumes isotropic materials. For birefringent media, additional checks may be needed.

## Tolerance Parameters

The internal `_validate_physics` function uses these default tolerances:

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `atol` | 1e-6 | Absolute tolerance for R + T ≈ 1 |
| `k_threshold` | 1e-10 | Threshold for considering a layer lossless |
