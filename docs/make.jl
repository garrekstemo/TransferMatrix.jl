push!(LOAD_PATH, "../src/")

using Documenter
using TransferMatrix

makedocs(
    sitename = "TransferMatrix.jl",
    modules = [TransferMatrix]
)