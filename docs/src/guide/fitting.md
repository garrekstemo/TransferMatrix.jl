# Fitting and Thick Substrates

Two features aimed at modelling measured spectra: automatic differentiation, for
gradient-based fits of thicknesses and indices, and incoherent (thick) layers, for
substrates and windows whose interference fringes the spectrometer does not resolve.

## Gradients with ForwardDiff

[`transfer`](@ref) and [`amplitude_coefficients`](@ref) are written generically in
their number type, so ForwardDiff dual numbers pass through layer thicknesses,
dispersion functions, Euler angles, permeabilities, sheet conductivities, `λ`, `θ`
and `μ`. The gradient of a transmittance with respect to the film thicknesses:

```julia
using TransferMatrix, ForwardDiff

λ = 1e4 / 2100                                   # 2100 cm⁻¹, in μm
stack(p) = [Layer(λ -> 1.0, 0.0),                # ambient
            Layer(λ -> 2.4, p[1]),               # film 1, thickness p[1]
            Layer(λ -> 1.35 + p[3]*im, p[2]),    # film 2, thickness p[2], extinction p[3]
            Layer(λ -> 1.45, 0.0)]               # substrate

ForwardDiff.gradient(p -> transfer(λ, stack(p); θ=0.3).Tpp, [0.5, 0.9, 0.01])
```

The test suite checks these gradients against central finite differences
(agreement better than 1e-6) on the 2×2 (`:auto`, `:decoupled`) and 4×4 (`:exp`)
paths, for isotropic, anisotropic, rotated and gyrotropic layers, sheets, the
circular basis and incoherent layers.

Not differentiable: `method=:eig` and an anisotropic or magnetic ambient or
substrate whose own parameters (or `θ`) carry dual numbers, because both
diagonalize with `eigen`. The sweeps and the field and energy functions return
`Float64`.

## Incoherent layers

A 3 mm CaF₂ window has etalon fringes about 1 cm⁻¹ apart, below the resolution of a
typical FTIR measurement, so the instrument records their average. Mark such layers
with `incoherent`:

```julia
air  = Layer(λ -> 1.0, 0.0)
caf2 = Layer(λ -> 1.40 + 1e-6im, 3000.0)         # 3 mm, in μm
H, L = Layer(λ -> 2.4, 0.5), Layer(λ -> 1.35, 0.9)

coated = [air, H, L, caf2, L, H, air]             # coatings on both faces
transfer(λ, coated; θ=0.2, incoherent=4)          # layer 4 is added by intensity
```

The stack is split at each incoherent layer into coherent blocks. Each block is
solved with the 2×2 path for its forward and backward power reflectance and
transmittance, and the blocks are chained through the incoherent layers'
single-pass attenuation with the intensity-matrix method of Katsidis and Siapkas
(Appl. Opt. 41, 3978, 2002). For a lossless thick layer the result equals the
coherent spectrum averaged over one period of the layer's round-trip phase; the
test suite checks this to 1e-12.

Several layers may be incoherent (`incoherent=[3, 7]`), and the keyword is also
accepted by [`sweep_angle`](@ref) and [`sweep_thickness`](@ref). Incoherent layers
must be isotropic and non-magnetic, every coherent block must decouple s and p, and
only `basis=:linear` is available; the cross-polarized terms are zero.
