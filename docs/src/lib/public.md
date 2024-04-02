# Public Documentation

These functions and types are to be used for transfer matrix calculation
based on the sources used. If you wish to modify any of the steps in the 
calculation, refer to the private API.

## Index

```@index
Pages = ["public.md"]
```

## Transfer Matrix Functions

```@autodocs
Modules = [TransferMatrix]
Pages = ["general_TMM.jl", "layer.jl", "matrix_constructors.jl"]
Private = false
```

## Miscellaneous Optics Functions

```@autodocs
Modules = [TransferMatrix]
Pages = ["optics_functions.jl"]
Private = false
```