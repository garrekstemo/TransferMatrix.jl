# Bloch surface wave (BSW) example
#
# A BSW is a surface mode at the termination of a truncated
# photonic crystal (DBR). It lives inside the photonic bandgap
# and appears as a sharp dip in the reflectance at a specific
# angle of incidence. The mode is excited via prism coupling
# (Kretschmann-like geometry): prism | DBR | air.

using RefractiveIndex
using TransferMatrix
using CairoMakie

# Materials
n_prism = RefractiveMaterial("glass", "soda-lime", "Rubin-clear")[1]
n_tio2 = RefractiveMaterial("main", "TiO2", "Sarkar")
n_sio2 = RefractiveMaterial("main", "SiO2", "Rodriguez-de_Marcos")
n_air = RefractiveMaterial("other", "air", "Ciddor")

# Quarter-wave stack at design wavelength
λ_0 = 0.633  # μm
t_tio2 = λ_0 / (4 * n_tio2(λ_0))
t_sio2 = λ_0 / (4 * n_sio2(λ_0))

# Truncated DBR — the top layer is TiO2 (high-index),
# which creates the surface termination that supports the BSW.
prism = Layer(n_prism, 0.5)
tio2 = Layer(n_tio2, t_tio2)
sio2 = Layer(n_sio2, t_sio2)
air = Layer(n_air, 0.5)

nperiods = 8
unit = [tio2, sio2]
layers = [prism, repeat(unit, nperiods)..., tio2, air]

# Angle-resolved reflectance
λs = range(0.50, 0.80, length=400)
θs = deg2rad.(range(40, 70, length=300))
spectra = sweep_angle(λs, θs, layers)

# Plot — the BSW appears as a narrow dark line
# threading through the bright stopband region.
fig = Figure(size=(600, 500))
ax = Axis(fig[1, 1],
    xlabel="Wavelength (nm)",
    ylabel="Angle of incidence (°)",
    title="Bloch Surface Wave in Photonic Bandgap",
)
hm = heatmap!(ax, λs .* 1e3, rad2deg.(θs), spectra.Rss',
    colormap=:inferno, colorrange=(0, 1))
Colorbar(fig[1, 2], hm, label="Rss")

outpath = joinpath(@__DIR__, "..", "docs", "src", "assets", "examples", "bloch_surface_wave.png")
mkpath(dirname(outpath))
save(outpath, fig)
println("saved ", outpath)
