# TransferMatrix.jl

[![codecov](https://codecov.io/gh/garrekstemo/TransferMatrix.jl/graph/badge.svg?token=WO2ITE125C)](https://codecov.io/gh/garrekstemo/TransferMatrix.jl)
[![CI](https://github.com/garrekstemo/TransferMatrix.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/garrekstemo/TransferMatrix.jl/actions/workflows/CI.yml)
[![stable](https://img.shields.io/badge/docs-stable-blue)](https://garrek.org/TransferMatrix.jl/stable/)

A general 4x4 transfer-matrix optics for layered media in Julia, with support for
isotropic stacks and polarization-resolved results.

## Highlights

- Isotropic multilayers with polarization-resolved reflection/transmission
- Angle and thickness sweeps for dispersion maps
- RefractiveIndex.jl integration plus custom dispersions
- Field profiles and interface positions for visualization
- Anisotropic layers and cross-polarization available but still experimental

## Installation

```julia
julia>]
pkg> add TransferMatrix
```

```julia
using TransferMatrix
```

## Quickstart

```julia
using TransferMatrix
using RefractiveIndex

λ = 0.6328  # μm
air = Layer(λ -> 1.0 + 0.0im, 0.0)
glass = Layer(RefractiveMaterial("main", "SiO2", "Malitson"), 0.5)

layers = [air, glass, air]
result = transfer(λ, layers; θ=0.0)
@show result.Rpp result.Tpp
```

## Anisotropic layer + angle sweep

```julia
using TransferMatrix

no = λ -> 1.658
ne = λ -> 1.486
crystal = Layer(no, no, ne, 0.5; euler=(0, pi/4, 0))
air = Layer(λ -> 1.0 + 0.0im, 0.0)

λs = range(0.45, 0.75, length=200)
θs = range(0.0, 0.8, length=150)
result = sweep_angle(λs, θs, [air, crystal, air])
Rpp = result.Rpp
```

## More examples

Check `examples/` for full scripts: DBR cavities, Brewster angle, Fabry-Perot,
polariton dispersion, and thickness sweeps.

## Documentation

For tutorials, API reference, and validation tests, see
https://garrek.org/TransferMatrix.jl.
