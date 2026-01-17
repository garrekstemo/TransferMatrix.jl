using Documenter, TransferMatrix

makedocs(
    sitename = "TransferMatrix.jl",
    modules = [TransferMatrix],
    pages = [
        "Introduction" => "index.md",
        "Guide" => Any[
                    "Quick Start" => "guide/quickstart.md",
                    "Tutorial" => "guide/tutorial.md",
                    "Physics Validation" => "guide/validation.md"
        ],
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