push!(LOAD_PATH, "../src/")

using Documenter, TransferMatrix

makedocs(
    sitename = "TransferMatrix.jl",
    modules = [TransferMatrix],
    pages = [
        "Introduction" => "index.md",
        "Tutorial" => "tutorial.md",
        "API Reference" => "reference.md",
        "References" => "bibliography.md"
    ]
)

deploydocs(
    repo = "github.com/garrekstemo/TransferMatrix.jl.git",
)