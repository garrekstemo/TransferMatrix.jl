push!(LOAD_PATH, "../src/")

using Documenter, TransferMatrix

makedocs(
    sitename = "TransferMatrix.jl",
    modules = [TransferMatrix],
)