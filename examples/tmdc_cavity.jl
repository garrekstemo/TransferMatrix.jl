# TMDC monolayer inside a Fabry–Pérot-like cavity, modeled as a conductive sheet.
# Demonstrates a Lorentzian-exciton sheet conductivity producing a polariton
# anticrossing in an angle sweep, plus an E/H field overlay.
#
# Run from the examples/ environment:
#   julia --project=examples examples/tmdc_cavity.jl

using TransferMatrix
using CairoMakie

# Lorentzian exciton sheet conductivity σ(λ) in SI Siemens.
#
# Under the package's exp(-iωt) convention an absorbing sheet has Re σ > 0, so the
# Lorentzian carries a leading −i:
#   σ(ω) = −i ε₀ c₀ f ω / (ω₀² − ω² − i ω Γ)
# On resonance (ω = ω₀) this gives σ_res = ε₀ c₀ f / Γ, which is real and positive.
#
# The peak strength is set in units of the universal sheet conductivity
# σ₀ = e²/(4ħ) ≈ 6.1e-5 S (dimensionless σ̃₀ = Z₀σ₀ ≈ 0.023), the natural scale for a
# TMDC monolayer. We choose the oscillator strength f so that σ_res = A·σ₀. A is set
# large here (A = 100) so the light–matter coupling clears the cavity/exciton
# linewidths and a genuine Rabi-split doublet appears.
function exciton_sigma(λ)
    c0 = 299792458.0
    ε0 = 8.8541878128e-12
    e  = 1.602176634e-19
    ħ  = 1.054571817e-34
    σ0 = e^2 / (4ħ)         # universal sheet conductivity scale (S)
    λ0 = 0.62               # exciton wavelength (μm)
    Γλ = 0.008              # linewidth (μm)
    # NOTE: the exciton coupling strength A here is an ILLUSTRATIVE knob chosen to give a
    # clearly visible splitting — NOT a measured WS2/MoSe2 oscillator strength. (And in
    # bsw_tmdc_polariton.jl the BSW crossing angle is hand-tuned.) These remain runnable
    # advanced examples; they are intentionally not part of the curated docs set.
    A  = 100.0              # peak conductivity in units of σ₀ (σ_res ≈ 100 σ₀)
    ω  = 2π * c0 / λ
    ω0 = 2π * c0 / λ0
    Γ  = 2π * c0 * Γλ / λ0^2
    f  = A * σ0 * Γ / (ε0 * c0)   # oscillator strength giving σ_res = A σ₀
    return -im * ε0 * c0 * f * ω / (ω0^2 - ω^2 - im * ω * Γ)
end

# Detuned half-wavelength cavity formed by two distributed Bragg mirrors. Each mirror
# is N quarter-wave (lo/hi) periods, giving enough finesse to resolve the polariton
# splitting. The high-index layer sits next to the spacer on both sides, so the
# half-wave spacer carries an electric-field antinode at its center — where the
# monolayer is placed. The spacer is tuned to λ_cav > λ₀ so the cavity mode starts red
# of the exciton at normal incidence and blueshifts through it as the angle increases,
# tracing an avoided crossing.
n_air, n_hi, n_lo, n_spacer = 1.0, 2.5, 1.45, 1.0
λ0     = 0.62                  # exciton wavelength (μm)
λ_cav  = 0.63                  # cavity-mode wavelength at normal incidence (μm)
N      = 6                     # quarter-wave periods per mirror

mirror_hi   = Layer(λ -> complex(n_hi), λ0 / (4 * n_hi))
mirror_lo   = Layer(λ -> complex(n_lo), λ0 / (4 * n_lo))
spacer_half = Layer(λ -> complex(n_spacer), λ_cav / (2 * n_spacer) / 2)
air = Layer(λ -> complex(n_air), 0.0)
sub = Layer(λ -> complex(n_air), 0.0)

# lo/hi periods so the high-index layer is adjacent to the spacer (field antinode there).
front = reduce(vcat, [[mirror_lo, mirror_hi] for _ in 1:N])
back  = reduce(vcat, [[mirror_hi, mirror_lo] for _ in 1:N])
layers = [air, front..., spacer_half, spacer_half, back..., sub]
# The two spacer halves are layers (1 + 2N + 1) and (1 + 2N + 2); their shared
# interface (the cavity-center antinode) is interface key 1 + 2N + 1.
sheet_index = 1 + 2N + 1                   # interface between the two spacer halves (cavity center)
sheets = Dict(sheet_index => Sheet(exciton_sigma))

λs = range(0.56, 0.69, length = 400)
θs = range(0.0, 0.5, length = 200)
spec = sweep_angle(collect(λs), collect(θs), layers; sheets = sheets)

fig = Figure(size = (900, 400))
ax1 = Axis(fig[1, 1], xlabel = "wavelength (μm)", ylabel = "angle (rad)", title = "Rss (polariton anticrossing)")
heatmap!(ax1, λs, θs, spec.Rss', colormap = :viridis)
vlines!(ax1, λ0, color = (:white, 0.6), linestyle = :dash)

# Field overlay computed for the BARE cavity (no sheet) at the exciton wavelength, to
# show that the monolayer is placed at a cavity-mode field antinode. With the sheet
# present the strong coupling opens the polariton gap and depletes the field exactly at
# λ₀, so the bare-cavity field is the one that makes the antinode placement visible.
ef = efield(λ0, layers; dz = 0.001)
hf = hfield(λ0, layers; dz = 0.001)
Emag = vec(sum(abs2, ef.s; dims = 1)) .^ 0.5
Hmag = vec(sum(abs2, hf.s; dims = 1)) .^ 0.5
# Interface key k sits at z = boundaries[k]; the cavity-center sheet is at boundaries[sheet_index].
z_sheet = ef.boundaries[sheet_index]       # z-position of the conductive sheet
ax2 = Axis(fig[1, 2], xlabel = "z (μm)", ylabel = "|field| (norm.)", title = "E / H at the monolayer")
lines!(ax2, ef.z, Emag ./ maximum(Emag), label = "|E| (s)")
lines!(ax2, hf.z, Hmag ./ maximum(Hmag), label = "|H| (s)")
for b in ef.boundaries
    vlines!(ax2, b, color = (:gray, 0.4))
end
vlines!(ax2, z_sheet, color = (:red, 0.8), linestyle = :dash, label = "sheet")
axislegend(ax2)

save(joinpath(@__DIR__, "tmdc_cavity.png"), fig)
println("wrote tmdc_cavity.png")
