# TMDC monolayer inside a Fabry–Pérot-like cavity, modeled as a conductive sheet.
# Demonstrates a Lorentzian-exciton sheet conductivity producing a polariton
# anticrossing in an angle sweep, plus an E/H field overlay.
#
# Run from the examples/ environment:
#   julia --project=examples examples/tmdc_cavity.jl

using TransferMatrix
using CairoMakie

# Lorentzian exciton sheet conductivity σ(λ) in SI Siemens.
# σ(ω) = i ε₀ c₀ f ω / (ω₀² - ω² - i ω Γ) with oscillator strength f (length units).
function exciton_sigma(λ)
    c0 = 299792458.0
    ε0 = 8.8541878128e-12
    λ0 = 0.62               # exciton wavelength (μm)
    Γλ = 0.01               # linewidth (μm)
    f = 5e-3                # oscillator strength (μm)
    ω  = 2π * c0 / λ
    ω0 = 2π * c0 / λ0
    Γ  = 2π * c0 * Γλ / λ0^2
    return im * ε0 * c0 * f * ω / (ω0^2 - ω^2 - im * ω * Γ)
end

# Simple half-wavelength cavity: two dielectric mirrors around a spacer, monolayer
# at the spacer center (a field antinode).
n_air, n_hi, n_spacer = 1.0, 2.2, 1.5
λ0 = 0.62
mirror_hi = Layer(λ -> complex(n_hi), λ0 / (4 * n_hi))
spacer_half = Layer(λ -> complex(n_spacer), λ0 / (2 * n_spacer) / 2)

air = Layer(λ -> complex(n_air), 0.0)
sub = Layer(λ -> complex(n_air), 0.0)
layers = [air, mirror_hi, spacer_half, spacer_half, mirror_hi, sub]
sheet_index = 3                          # interface between the two spacer halves (cavity center)
sheets = Dict(sheet_index => Sheet(exciton_sigma))

λs = range(0.58, 0.66, length = 400)
θs = range(0.0, 0.6, length = 200)
spec = sweep_angle(collect(λs), collect(θs), layers; sheets = sheets)

fig = Figure(size = (900, 400))
ax1 = Axis(fig[1, 1], xlabel = "wavelength (μm)", ylabel = "angle (rad)", title = "Rss (polariton anticrossing)")
heatmap!(ax1, λs, θs, spec.Rss', colormap = :viridis)

ef = efield(λ0, layers; dz = 0.001, sheets = sheets)
hf = hfield(λ0, layers; dz = 0.001, sheets = sheets)
ax2 = Axis(fig[1, 2], xlabel = "z (μm)", ylabel = "|field| (norm.)", title = "E / H at the monolayer")
lines!(ax2, ef.z, vec(sum(abs2, ef.s; dims = 1)) .^ 0.5, label = "|E| (s)")
lines!(ax2, hf.z, vec(sum(abs2, hf.s; dims = 1)) .^ 0.5, label = "|H| (s)")
for b in ef.boundaries
    vlines!(ax2, b, color = (:gray, 0.4))
end
axislegend(ax2)

save("tmdc_cavity.png", fig)
println("wrote tmdc_cavity.png")
