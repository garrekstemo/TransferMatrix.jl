# Regenerate every committed docs figure by running the curated examples.
# Each curated example saves its PNG to docs/src/assets/examples/<name>.png.
#   julia --project=examples examples/render_figures.jl

curated = [
    "fabry-perot",
    "dbr_cavity",
    "bloch_surface_wave",
    "graphene_absorption",
    "polariton_dispersion",
    "cholesteric_circular_bragg",
    "anisotropic_energy_budget",
    "thickness_dependence",
    "field_profiles",
]

for name in curated
    println("\n========== rendering ", name, " ==========")
    include(joinpath(@__DIR__, name * ".jl"))
end
println("\nAll ", length(curated), " curated figures rendered.")
