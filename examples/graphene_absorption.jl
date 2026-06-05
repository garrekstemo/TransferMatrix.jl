# Graphene's famous 2.3% — a one-atom-thick sanity check for conductive sheets.
#
# A single free-standing graphene sheet absorbs π·α ≈ 2.3% of visible light — a value
# fixed by the fine-structure constant α alone, and essentially the same at every
# visible color (R. R. Nair et al., "Fine structure constant defines visual
# transparency of graphene", Science 320, 1308 (2008)).
#
# Here we model graphene by its universal optical conductivity σ₀ = e²/4ħ (the real,
# absorptive part, which is accurate across the visible) as a zero-thickness `Sheet`,
# and check that TransferMatrix.jl reproduces that ~2.3% absorption. It is about the
# simplest possible test that the conductive-sheet feature does the right thing: the
# answer is a famous number you can look up.
#
# Run from the examples environment. The examples project normally uses the registered
# TransferMatrix, so point it at this checkout first:
#   julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia --project=examples examples/graphene_absorption.jl

using TransferMatrix
using CairoMakie

# --- physical constants ---
e = 1.602176634e-19       # elementary charge (C)
ħ = 1.054571817e-34       # reduced Planck constant (J·s)
α = 1 / 137.035999084     # fine-structure constant (dimensionless)
σ₀ = e^2 / (4ħ)           # graphene's universal optical conductivity (siemens)

# A "free-standing" sheet just means vacuum on both sides. The two bounding media are
# semi-infinite, so their thickness (0.0 here) is ignored by the calculation.
vacuum = Layer(λ -> 1.0 + 0.0im, 0.0)
layers = [vacuum, vacuum]

# Place the graphene at the one interface (between layer 1 and layer 2). A real,
# positive conductivity makes the sheet absorptive.
sheets = Dict(1 => Sheet(σ₀))

# Absorbed fraction at wavelength λ (μm): whatever is neither transmitted nor reflected.
function absorption(λ)
    r = transfer(λ, layers; sheets = sheets)
    return 1 - r.Tss - r.Rss
end

# --- the check, at green light (550 nm) ---
λ = 0.55
r = transfer(λ, layers; sheets = sheets)
T, R = r.Tss, r.Rss
A = 1 - T - R

println("One graphene sheet at λ = 550 nm:")
println("  transmitted  T = $(round(100T, digits = 2)) %")
println("  reflected    R = $(round(100R, digits = 3)) %")
println("  absorbed     A = $(round(100A, digits = 2)) %")
println("  textbook  π·α  = $(round(100π * α, digits = 2)) %   (Nair et al., Science 2008)")
println("  energy R+T+A   = $(round(100 * (R + T + A), digits = 4)) %   (conserved → 100)")
println(isapprox(A, π * α; rtol = 0.05) ?
        "  ✓ reproduces graphene's famous ~2.3% absorption" :
        "  ✗ unexpected: does not match π·α")

# --- and it holds across the whole visible spectrum ---
λs = range(0.40, 0.78, length = 200)      # μm, violet → red
Aλ = absorption.(λs)

fig = Figure(size = (680, 420))
ax = Axis(fig[1, 1];
    xlabel = "wavelength (μm)",
    ylabel = "absorption (%)",
    title = "A single graphene sheet absorbs ≈ 2.3% at every visible color")
lines!(ax, λs, 100 .* Aλ; linewidth = 3, label = "TransferMatrix.jl")
hlines!(ax, 100π * α; color = :red, linestyle = :dash, label = "π·α (textbook, Nair 2008)")
ylims!(ax, 0, 3)
axislegend(ax; position = :rb)
save("graphene_absorption.png", fig)
println("wrote graphene_absorption.png")
