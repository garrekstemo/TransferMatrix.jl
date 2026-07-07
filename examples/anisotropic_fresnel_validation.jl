# Oracle validation of the per-eigenmode cross-polarized transmittance.
#
# transfer() reports four transmittances per input polarization: the
# co-polarized Tpp, Tss and the converted Tps, Tsp. Each is the Poynting flux
# of ONE substrate eigenmode — evaluated with that mode's own wavevector —
# divided by the incident flux. The energy-budget example
# (anisotropic_energy_budget.jl) shows that these four channels close
# Rpp + Rps + Tpp + Tps = 1 internally; a budget can close, however, and still
# split the power wrongly between the channels. Here each channel is checked
# against independent closed-form textbook physics the Berreman machinery was
# never fitted to.
#
# Validation 1 — per-axis Fresnel oracle (quantitative, < 1e-10).
#   Air over a semi-infinite uniaxial crystal with the optic axis along ŷ:
#   Layer(n_o, n_e, n_o, ...). With the plane of incidence x–z, the two
#   transmitted eigenmodes are PURE p and s waves: the p mode (E in x–z) sees
#   only ε_xx = ε_zz = n_o², the s mode (E ∥ ŷ) only ε_yy = n_e². So Tpp must
#   equal the isotropic single-interface p-polarization Fresnel transmittance
#   with index n_o, and Tss the s-polarization one with index n_e — including
#   the flux factor n·cosθ_t/cosθ_i, which DIFFERS between the two output
#   modes because each refracts with its own index. Cross terms must vanish.
#
# Validation 2 — half-wave plate oracle (quantitative, < 1e-10).
#   A calcite-like slab (n_o = 1.658, n_e = 1.486) with the optic axis in the
#   layer plane at 45° to the plane of incidence (euler = (π/4, π/2, 0)) and
#   thickness t = λ₀/(2Δn) at λ₀ = 0.59 μm: a half-wave plate. At normal
#   incidence a p input splits equally onto the two eigenaxes; the slab
#   retards one by π and the recombined polarization is rotated by 90°, so the
#   transmitted power exits almost entirely s-polarized (Tps). Two formulas:
#     (a) the classic waveplate conversion T_conv = sin²(2·45°)·sin²(πΔn·t/λ),
#         scaled by the single-pass Fresnel surface losses (fringe-free);
#     (b) an EXACT Jones × Airy oracle — at normal incidence each eigenmode is
#         an independent Fabry–Pérot etalon with its own index, so
#         Tps = |t_e − t_o|²/4 and Tpp = |t_e + t_o|²/4 in closed form.
#   (b) must agree with transfer() to < 1e-10 across the sweep; (a) is the
#   envelope the etalon fringes oscillate around.
#
# Validation 3 — the reading the fix changed (for readers migrating).
#   In releases ≤ v3.2.1, Tpp was the TOTAL per-input transmitted flux and
#   Tps = |t_ps|² a raw amplitude ratio: at the half-wave point
#   Rpp + Rps + Tpp = 1 looked like a closed budget (with Tpp ≈ 0.848), while
#   the four-channel sum double-counted the converted power to ≈ 1.85. Now
#   Tpp is the co-polarized SHARE (≈ 0.0004 here) and the budget that closes
#   is the four-channel one.
#
# Run from the examples environment. The examples project normally uses the
# registered TransferMatrix, so point it at this checkout first:
#   julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia --project=examples examples/anisotropic_fresnel_validation.jl

using TransferMatrix
using CairoMakie

function check_below(label, deviation, tol)
    ok = deviation < tol
    println("  ", label, " = ", deviation, ok ? "  (PASS < $tol)" : "  (FAIL, tol $tol)")
    ok || error("validation failed: " * label)
end

no, ne = 1.658, 1.486
λ0 = 0.59
air = Layer(λ -> 1.0, 1.0)

# --- Validation 1: air | uniaxial half-space, optic axis along ŷ ------------
# Textbook single-interface Fresnel transmittance (incident index n₁ = 1),
# written out for a substrate of index n:
#   cosθ_t = √(1 − sin²θᵢ/n²)                       (Snell)
#   t_p = 2cosθᵢ/(n·cosθᵢ + cosθ_t)                 (amplitude)
#   t_s = 2cosθᵢ/(cosθᵢ + n·cosθ_t)
#   T   = (n·cosθ_t/cosθᵢ)·|t|²                     (transmitted flux ratio)
# The p mode refracts with n_o and the s mode with n_e, so each channel
# carries its OWN n·cosθ_t flux factor.

substrate_y = Layer(λ -> no, λ -> ne, λ -> no, 1.0)
θdeg = 0.0:1.0:80.0
θi = deg2rad.(θdeg)
res1 = [transfer(λ0, [air, substrate_y]; θ = θ) for θ in θi]
Tpp1 = getfield.(res1, :Tpp)
Tss1 = getfield.(res1, :Tss)

cosθp = sqrt.(1 .- (sin.(θi) ./ no) .^ 2)
cosθs = sqrt.(1 .- (sin.(θi) ./ ne) .^ 2)
tp = 2 .* cos.(θi) ./ (no .* cos.(θi) .+ cosθp)
ts = 2 .* cos.(θi) ./ (cos.(θi) .+ ne .* cosθs)
Tp_fresnel = no .* cosθp ./ cos.(θi) .* abs2.(tp)
Ts_fresnel = ne .* cosθs ./ cos.(θi) .* abs2.(ts)

cross1 = maximum(res -> max(res.Tps, res.Tsp, res.Rps, res.Rsp), res1)

println("Validation 1: air | uniaxial half-space (optic axis ∥ ŷ), θ = 0–80°")
check_below("max |Tpp − Fresnel T_p(n_o)|", maximum(abs.(Tpp1 .- Tp_fresnel)), 1e-10)
check_below("max |Tss − Fresnel T_s(n_e)|", maximum(abs.(Tss1 .- Ts_fresnel)), 1e-10)
check_below("max cross channel (Tps, Tsp, Rps, Rsp)", cross1, 1e-10)

# --- Validation 2: calcite-like half-wave plate at 45° ----------------------
Δn = no - ne
t_hwp = λ0 / (2Δn)
hwp = Layer(λ -> no, λ -> no, λ -> ne, t_hwp; euler = (π/4, π/2, 0.0))
λs = 0.40:0.002:1.00
res2 = [transfer(λ, [air, hwp, air]) for λ in λs]
Tps2 = getfield.(res2, :Tps)
Tpp2 = getfield.(res2, :Tpp)

# (a) Classic waveplate conversion, scaled by surface losses.
# a_j = t_in·t_out = [2/(1+n_j)]·[2n_j/(1+n_j)] is the two-face single-pass
# amplitude transmission of eigenmode j; the fringe-free envelope of the
# converted power is a_o·a_e·sin²(2β)·sin²(πΔn·t/λ) with β = 45°.
a_o = 4no / (1 + no)^2
a_e = 4ne / (1 + ne)^2
Tps_classic = a_o * a_e * sin(2 * (π/4))^2 .* sin.(π * Δn * t_hwp ./ λs) .^ 2

# (b) Exact Jones × Airy oracle. Each eigenmode traverses the slab as its own
# Fabry–Pérot etalon:
#   t_slab(n) = t_in·t_out·e^{iδ}/(1 − r²e^{2iδ}),  δ = 2πnt/λ,  r = (n−1)/(n+1)
# With the eigenaxes at 45°, an x̂ (p) input recombines in the output air as
# E_x = (t_e + t_o)/2 and E_y = (t_e − t_o)/2, so
#   Tpp = |t_e + t_o|²/4,  Tps = |t_e − t_o|²/4.
function slab_transmission(n, t, λ)
    δ = 2π * n * t / λ
    r = (n - 1) / (n + 1)
    return (2 / (1 + n)) * (2n / (1 + n)) * cis(δ) / (1 - r^2 * cis(2δ))
end
t_o = slab_transmission.(no, t_hwp, λs)
t_e = slab_transmission.(ne, t_hwp, λs)
Tps_airy = abs2.(t_e .- t_o) ./ 4
Tpp_airy = abs2.(t_e .+ t_o) ./ 4

budget_p = [r.Rpp + r.Rps + r.Tpp + r.Tps for r in res2]
budget_s = [r.Rss + r.Rsp + r.Tss + r.Tsp for r in res2]

println("Validation 2: half-wave plate (t = ", round(t_hwp, digits = 5),
        " μm), λ = 0.40–1.00 μm, normal incidence")
check_below("max |Tps − Jones×Airy oracle|", maximum(abs.(Tps2 .- Tps_airy)), 1e-10)
check_below("max |Tpp − Jones×Airy oracle|", maximum(abs.(Tpp2 .- Tpp_airy)), 1e-10)
check_below("max |p-input budget − 1|", maximum(abs.(budget_p .- 1)), 1e-12)
check_below("max |s-input budget − 1|", maximum(abs.(budget_s .- 1)), 1e-12)

# The fringe-free conversion envelope peaks exactly at the design wavelength;
# the exact curve's global maximum rides the etalon fringe nearest to it
# (fringe spacing λ²/(2n̄t) ≈ 0.065 μm here).
i0 = argmin(abs.(λs .- λ0))
λpeak_env = λs[argmax(Tps_classic)]
λpeak = λs[argmax(Tps2)]
println("  conversion envelope peaks at λ = ", λpeak_env, " μm  (design λ₀ = ", λ0, " μm)")
println("  exact Tps max at λ = ", λpeak, " μm (nearest etalon fringe);  Tps(λ₀) = ",
        round(Tps2[i0], digits = 6), ",  Tpp(λ₀) = ", round(Tpp2[i0], digits = 6))
abs(λpeak_env - λ0) < 1e-9 || error("validation failed: envelope does not peak at the design wavelength")
abs(λpeak - λ0) < 0.05 || error("validation failed: Tps peak is not within one fringe of the design wavelength")
Tps2[i0] > 0.8 || error("validation failed: conversion at the design wavelength should be nearly complete")
Tpp2[i0] < 1e-3 || error("validation failed: Tpp at the half-wave point should be ≈ 0")

# --- Validation 3: the reading the per-mode fix changed ---------------------
# Migration note (releases ≤ v3.2.1 → this version): Tpp used to be the TOTAL
# transmitted Poynting flux for a p input (both output polarizations, one
# shared wavevector) and Tps = |t_ps|² a raw amplitude ratio. Rpp + Rps + Tpp
# was therefore the sum that read 1.0, and adding Tps double-counted the
# converted power (≈ 1.85 for this stack). With the per-eigenmode
# decomposition, Tpp at the half-wave point drops from ≈ 0.848 (the total) to
# ≈ 0.0004 (the co-polarized share), Rpp + Rps + Tpp is correctly ≪ 1, and
# the budget that closes is the four-channel one.
res_hw = transfer(λ0, [air, hwp, air])
old_reading = res_hw.Rpp + res_hw.Rps + res_hw.Tpp
println("Validation 3: half-wave point, per-mode decomposition vs the pre-fix reading")
println("  Tpp = ", round(res_hw.Tpp, digits = 6),
        "   (co-polarized share; ≤ v3.2.1 put the TOTAL ≈ 0.848 here)")
println("  Tps = ", round(res_hw.Tps, digits = 6),
        "   (converted share; ≤ v3.2.1 reported |t_ps|² instead)")
println("  Rpp + Rps + Tpp       = ", round(old_reading, digits = 6),
        "   (the old \"closed\" budget — correctly ≠ 1 now)")
println("  Rpp + Rps + Tpp + Tps = ", res_hw.Rpp + res_hw.Rps + res_hw.Tpp + res_hw.Tps,
        "   (the budget that closes)")
old_reading < 0.9 || error("validation failed: Tpp looks like a total flux again (pre-fix semantics?)")

# --- Figure ------------------------------------------------------------------
fig = Figure(size = (980, 430))

axa = Axis(fig[1, 1];
    xlabel = "Incidence angle (°)", ylabel = "Transmittance",
    title = "Air | uniaxial half-space: per-axis Fresnel oracle")
lines!(axa, θdeg, Tp_fresnel; color = :black, linewidth = 2.5,
    label = "Fresnel T_p(n_o) analytic")
lines!(axa, θdeg, Ts_fresnel; color = :gray55, linewidth = 2.5,
    label = "Fresnel T_s(n_e) analytic")
scatter!(axa, θdeg[1:4:end], Tpp1[1:4:end]; color = :orangered3, markersize = 8,
    label = "transfer() Tpp")
scatter!(axa, θdeg[1:4:end], Tss1[1:4:end]; color = :dodgerblue3, marker = :utriangle,
    markersize = 8, label = "transfer() Tss")
axislegend(axa; position = :lb, framevisible = true)

axb = Axis(fig[1, 2];
    xlabel = "Wavelength (μm)", ylabel = "Transmittance (p input)",
    title = "Half-wave plate at 45°: waveplate oracle")
limits!(axb, 0.38, 1.02, -0.04, 1.38)
lines!(axb, λs, Tps_classic; color = :gray55, linewidth = 2, linestyle = :dash,
    label = "classic sin²(πΔn·t/λ) × surface loss")
lines!(axb, λs, Tps2; color = :firebrick3, linewidth = 2.5, label = "transfer() Tps")
lines!(axb, λs, Tpp2; color = :dodgerblue3, linewidth = 1.8, label = "transfer() Tpp")
scatter!(axb, λs[1:15:end], Tps_airy[1:15:end]; color = :black, markersize = 6,
    label = "Jones × Airy exact")
vlines!(axb, λ0; color = :black, linestyle = :dot)
text!(axb, λ0 + 0.015, 0.06; text = "design λ₀", align = (:left, :center), fontsize = 11)
axislegend(axb; position = :rt, framevisible = true)

save(joinpath(@__DIR__, "anisotropic_fresnel_validation.png"), fig)
println("saved ", joinpath(@__DIR__, "anisotropic_fresnel_validation.png"))
