# Bloch surface wave coupled to a TMDC exciton — a BSW-exciton-polariton anticrossing.
#
# A Bloch surface wave (BSW) is a guided mode that lives at the truncated surface of a
# distributed Bragg reflector (DBR). Because its frequency sits inside the photonic
# bandgap of the DBR, the mode cannot radiate into the stack; it is bound to the
# surface and decays into the air above. To excite it we shine light through a glass
# prism beyond the critical angle (a Kretschmann/attenuated-total-reflection geometry):
#
#       prism  |  DBR (8 × TiO2/SiO2 periods + TiO2 cap)  |  air
#
# Beyond the critical angle no light is transmitted, so all-but-the-BSW is reflected and
# the BSW shows up as a sharp DIP in the reflectance Rss at the angle/wavelength where
# the surface mode is phase-matched. (Reflectance R = |r|² is exact in this total-
# internal-reflection regime — this is just ATR — so we never need transmittance here.)
#
# Now place a transition-metal-dichalcogenide (TMDC) monolayer (e.g. WS2/MoSe2) right at
# the DBR/air surface, exactly where the BSW field is largest. A monolayer is one atom-
# layer thick, so we do NOT model it as a bulk Layer; instead we model it as a 2D surface
# conductivity σ(λ) — the package's `Sheet` — carrying a Lorentzian exciton resonance.
# When the bare BSW dispersion is tuned onto the exciton, the photon-like BSW and the
# material exciton hybridize: the single BSW dip splits into TWO branches (an upper and a
# lower BSW-exciton-polariton) that bend away from the exciton wavelength λ₀ and never
# cross it. That avoided crossing (anticrossing) is the signature of strong light–matter
# coupling, and its size is the Rabi splitting.
#
# Run from the examples environment. The examples project normally uses the registered
# TransferMatrix, so point it at this checkout first:
#   julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia --project=examples examples/bsw_tmdc_polariton.jl

using RefractiveIndex
using TransferMatrix
using CairoMakie

# --- Materials -------------------------------------------------------------------------
# Glass prism (high index, to reach incidence angles beyond the DBR/air critical angle),
# the two DBR materials, and air on top.
n_prism = RefractiveMaterial("glass", "soda-lime", "Rubin-clear")[1]
n_tio2  = RefractiveMaterial("main", "TiO2", "Sarkar")
n_sio2  = RefractiveMaterial("main", "SiO2", "Rodriguez-de_Marcos")
n_air   = RefractiveMaterial("other", "air", "Ciddor")

# --- DBR geometry (quarter-wave stack at the design wavelength) ------------------------
λ_design = 0.633                              # μm, quarter-wave design wavelength
t_tio2 = λ_design / (4 * n_tio2(λ_design))    # high-index quarter-wave thickness (μm)
t_sio2 = λ_design / (4 * n_sio2(λ_design))    # low-index  quarter-wave thickness (μm)

prism = Layer(n_prism, 0.5)                   # semi-infinite incident medium
tio2  = Layer(n_tio2, t_tio2)
sio2  = Layer(n_sio2, t_sio2)
cap   = Layer(n_tio2, 0.5 * t_tio2)           # truncation cap: a half-thickness TiO2
air   = Layer(n_air, 0.5)                      # semi-infinite ambient

# 8 TiO2/SiO2 periods, then the thin TiO2 cap. The cap truncates the crystal so a surface
# state (the BSW) appears in the bandgap; its field peaks at the cap/air interface.
layers = [prism, repeat([tio2, sio2], 8)..., cap, air]   # 19 layers

# The conductive sheet goes at the cap/air surface — the interface just before air, which
# is where the BSW field is maximal (best overlap with the monolayer).
sheet_index = lastindex(layers) - 1            # = 18, the DBR/air (BSW field-max) surface

# --- TMDC exciton as a Lorentzian surface conductivity --------------------------------
# σ(λ) in SI Siemens. Under the package's exp(-iωt) convention an absorbing sheet has
# Re σ > 0, so the Lorentzian carries a leading −i:
#   σ(ω) = −i ε₀ c₀ f ω / (ω₀² − ω² − i ω Γ)
# On resonance σ_res = ε₀ c₀ f / Γ is real and positive. We scale the peak to the
# universal monolayer conductivity σ₀ = e²/(4ħ) ≈ 6.1e-5 S by choosing the oscillator
# strength f so that σ_res ≈ A·σ₀. A is the coupling knob: larger A → larger splitting.
function exciton_sigma(λ; λ0 = 0.5866, Γλ = 0.006, A = 30.0)
    c0 = 299792458.0
    ε0 = 8.8541878128e-12
    e  = 1.602176634e-19
    ħ  = 1.054571817e-34
    σ0 = e^2 / (4ħ)                  # universal sheet conductivity scale (S)
    ω  = 2π * c0 / λ
    ω0 = 2π * c0 / λ0
    Γ  = 2π * c0 * Γλ / λ0^2         # angular-frequency linewidth from Γλ (μm)
    f  = A * σ0 * Γ / (ε0 * c0)      # oscillator strength giving σ_res ≈ A·σ₀
    return -im * ε0 * c0 * f * ω / (ω0^2 - ω^2 - im * ω * Γ)
end

# Exciton tuned onto the bare BSW dispersion; coupling A = 30 (peak σ in units of σ₀).
λ0_exc  = 0.5866                                # exciton wavelength (μm)
# NOTE: the exciton coupling strength A here is an ILLUSTRATIVE knob chosen to give a
# clearly visible splitting — NOT a measured WS2/MoSe2 oscillator strength. (And in
# bsw_tmdc_polariton.jl the BSW crossing angle is hand-tuned.) These remain runnable
# advanced examples; they are intentionally not part of the curated docs set.
A_final = 30.0                                  # coupling strength (σ_res ≈ 30 σ₀)
sheets  = Dict(sheet_index => Sheet(λ -> exciton_sigma(λ; λ0 = λ0_exc, A = A_final)))

# --- Angle-resolved reflectance maps (s-polarization) ---------------------------------
# Same windows for both maps. Rss dips trace the BSW (bare) and the polaritons (coupled).
λ_window = range(0.566, 0.606, length = 1001)              # μm
θ_window = range(deg2rad(44.5), deg2rad(48.5), length = 161)  # radians

bare    = sweep_angle(collect(λ_window), collect(θ_window), layers)
coupled = sweep_angle(collect(λ_window), collect(θ_window), layers; sheets = sheets)
# bare.Rss    : a single BSW dark line crossing the exciton wavelength λ₀.
# coupled.Rss : the line splits into upper/lower polariton branches that avoid λ₀.

# --- Side-by-side heatmaps ------------------------------------------------------------
# Rss has shape (length(θ), length(λ)); heatmap(x=λ, y=θ) wants z as (length(λ), length(θ)),
# so we transpose. Shared color scale 0..1 makes the two panels directly comparable.
λ_nm    = λ_window .* 1e3                       # nm for the x-axis
θ_deg   = rad2deg.(θ_window)                    # degrees for the y-axis
λ0_nm   = λ0_exc * 1e3

fig = Figure(size = (1000, 440))

ax1 = Axis(fig[1, 1],
    xlabel = "wavelength (nm)",
    ylabel = "angle of incidence (°)",
    title  = "Bare BSW")
heatmap!(ax1, λ_nm, θ_deg, bare.Rss', colormap = :inferno, colorrange = (0, 1))
vlines!(ax1, λ0_nm, color = (:cyan, 0.8), linestyle = :dash)

ax2 = Axis(fig[1, 2],
    xlabel = "wavelength (nm)",
    ylabel = "angle of incidence (°)",
    title  = "BSW + TMDC monolayer (anticrossing)")
hm = heatmap!(ax2, λ_nm, θ_deg, coupled.Rss', colormap = :inferno, colorrange = (0, 1))
vlines!(ax2, λ0_nm, color = (:cyan, 0.8), linestyle = :dash)

Colorbar(fig[1, 3], hm, label = "Rss")

save("bsw_tmdc_polariton.png", fig)
println("wrote bsw_tmdc_polariton.png")

# --- Cleanest single-trace proof: a wavelength cut at the bare crossing angle ----------
# At the angle where the bare BSW crosses λ₀, the bare reflectance has ONE dip at λ₀,
# while the coupled reflectance has TWO dips straddling λ₀ — the polariton doublet.
θ_cross = deg2rad(46.38)
λ_cut   = range(0.560, 0.610, length = 2501)
R_bare  = sweep_angle(collect(λ_cut), [θ_cross], layers).Rss[1, :]
R_cpl   = sweep_angle(collect(λ_cut), [θ_cross], layers; sheets = sheets).Rss[1, :]

# Find reflectance dips in a 1D trace by prominence: a local minimum counts if it sits at
# least `prom` below the lower of the two surrounding local maxima. This catches the
# shallow lower-polariton dip as well as the deep upper one (a fixed global-range cutoff
# would miss the shallow branch).
function find_dips(R; prom = 0.02)
    minima = [i for i in 2:(length(R) - 1) if R[i] < R[i - 1] && R[i] < R[i + 1]]
    dips = Int[]
    for i in minima
        left  = maximum(@view R[1:i])          # highest point on the way in
        right = maximum(@view R[i:end])         # highest point on the way out
        if min(left, right) - R[i] ≥ prom
            push!(dips, i)
        end
    end
    return dips
end

bare_dips = find_dips(R_bare)
cpl_dips  = find_dips(R_cpl)
λ_bare    = collect(λ_cut)[bare_dips]
λ_cpl     = collect(λ_cut)[cpl_dips]

println()
println("Wavelength cut at θ = $(round(rad2deg(θ_cross), digits = 2))° (bare BSW crossing angle):")
println("  exciton wavelength λ₀ = $(round(λ0_exc * 1e3, digits = 1)) nm,  coupling A = $(A_final) σ₀")
println("  bare BSW dips   : ", round.(λ_bare .* 1e3, digits = 1), " nm   (", length(bare_dips), " dip)")
println("  coupled dips    : ", round.(λ_cpl .* 1e3, digits = 1), " nm   (", length(cpl_dips), " dips)")
if length(cpl_dips) == 2
    splitting = abs(λ_cpl[2] - λ_cpl[1])
    println("  → Rabi splitting ≈ $(round(splitting * 1e3, digits = 1)) nm  ($(round(splitting, digits = 4)) μm)")
    println("  ✓ anticrossing visible: one bare BSW dip splits into two polariton branches")
else
    println("  ✗ unexpected: did not resolve a two-dip polariton doublet")
end
