# Internals

These functions are not exported by the TransferMatrix.jl module
and can be called using the `TransferMatrix.` qualifier. Use
these methods if you wish to construct a transfer matrix method
manually step by step or modify intermediate steps.

## Index

```@index
Pages = ["internals.md"]
```

## General Transfer Matrix Method

```@autodocs
Modules = [TransferMatrix]
Pages = ["general_TMM.jl", "layer.jl", "matrix_constructors.jl"]
Public = false
```