push!(LOAD_PATH, "../src/")

using Documenter, TransferMatrix

makedocs(
    sitename = "TransferMatrix.jl",
    modules = [TransferMatrix],
    pages = [
        "Introduction" => "index.md",
        "Tutorial" => "tutorial.md",
        "Library" => Any[
                        "Public" => "lib/public.md",
                        "Internals" => "lib/internals.md"
        ],
        "References" => "bibliography.md"
    ]
)

deploydocs(
    repo = "github.com/garrekstemo/TransferMatrix.jl.git",
)