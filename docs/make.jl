using Documenter, TransferMatrix

makedocs(
    sitename = "TransferMatrix.jl",
    modules = [TransferMatrix],
    pages = [
        "Introduction" => "index.md",
        "Guide" => Any[
                    "Quick Start" => "guide/quickstart.md",
                    "Tutorial" => "guide/tutorial.md",
                    "Dispersion Models" => "guide/dispersion_models.md",
                    "Physics Validation" => "guide/validation.md"
        ],
        "Examples" => Any[
                    "Fabry–Pérot cavity" => "examples/fabry_perot.md",
                    "DBR microcavity" => "examples/dbr_cavity.md",
                    "Bloch surface wave" => "examples/bloch_surface_wave.md",
                    "Graphene absorption" => "examples/graphene_absorption.md",
                    "Vibrational polariton dispersion" => "examples/polariton_dispersion.md",
                    "Cholesteric circular Bragg" => "examples/cholesteric_circular_bragg.md",
                    "Anisotropic energy budget" => "examples/anisotropic_energy_budget.md",
                    "Thickness dependence" => "examples/thickness_dependence.md",
                    "Field profiles (E and H)" => "examples/field_profiles.md",
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