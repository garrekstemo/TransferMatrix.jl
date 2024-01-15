# Internals

These functions are not exported by the TransferMatrix.jl module
and can be called using the `TransferMatrix.` qualifier. Use
these methods if you wish to construct a transfer matrix method
manually step by step or modify intermediate steps.

## Index

```@index
Pages = ["internals.md"]
```

## Transfer Matrix Functions

```@autodocs
Modules = [TransferMatrix]
Pages = ["general_TMM.jl", "layers.jl"]
Public = false
```

## Types

```@autodocs
Modules = [TransferMatrix]
Pages = ["types.jl"]
Public = false
```

## Data Read/Write Functions

```@autodocs
Modules = [TransferMatrix]
Pages = ["dataio.jl"]
Public = false
```

## Optics Functions

```@autodocs
Modules = [TransferMatrix]
Pages = ["optics_functions.jl"]
Public = false
```