# TransferMatrix.jl

[![codecov](https://codecov.io/gh/garrekstemo/TransferMatrix.jl/graph/badge.svg?token=WO2ITE125C)](https://codecov.io/gh/garrekstemo/TransferMatrix.jl)
[![CI](https://github.com/garrekstemo/TransferMatrix.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/garrekstemo/TransferMatrix.jl/actions/workflows/CI.yml)
[![stable](https://img.shields.io/badge/docs-stable-blue)](https://garrek.org/TransferMatrix.jl/stable/)

A general 4x4 transfer matrix method implementation for optical waves in layered media in Julia.

## Installation

```
julia>]
pkg> add TransferMatrix
```

Start using TransferMatrix.jl

```
using TransferMatrix
```

## Documentation

For more details, including a comprehensive tutorial, see the [documentation website](https://garrek.org/TransferMatrix.jl).


## TODOs

- [] Improve test coverage
- [] Use [dispersion](https://stillyslalom.github.io/RefractiveIndex.jl/dev/#RefractiveIndex.dispersion-Tuple%7BRefractiveMaterial,%20Float64%7D) to calculate refractive indices.
- [] `Layer` struct has only one refractive index (instead of an array for several wavelength), then propagate light of one wavelength at a time through the structure.
- [] Clean up types in function arguments.