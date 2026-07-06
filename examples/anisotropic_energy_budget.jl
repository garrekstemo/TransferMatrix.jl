# Energy conservation in a rotated anisotropic slab: the missing cross-polarization.
#
# A lossless dielectric slab must obey R + T = 1. But take a uniaxial crystal and
# rotate its optic axis OUT of the plane of incidence, and the naive budget
#
#     R + T  =  Rpp + Tpp  ≈  0.987
#
# falls short of 1 — even though nothing absorbs. The slab is not leaking energy:
# a rotated anisotropic crystal CONVERTS polarization. A purely p-polarized input
# comes back with a small s-polarized component in BOTH directions — the
# cross-polarized reflectance Rps and the cross-polarized transmittance Tps.
# Those converted channels are exactly what the naive budget forgets. The
# complete, polarization-resolved budget
#
#     Rpp + Rps + Tpp + Tps  =  1
#
# closes to machine precision. (Each transmittance is the Poynting flux of its
# own substrate eigenmode, so the four channels partition the energy exactly.)
#
# Panel (a): the naive vs. complete budget vs. incidence angle — the shaded gap
#            IS the converted power Rps + Tps.
# Panel (b): the converted channels vs. the optic-axis azimuth — conversion
#            vanishes when the axis lies in the plane of incidence and grows as
#            it tilts out.
#
# Run from the examples environment:
#   julia --project=examples examples/anisotropic_energy_budget.jl

using TransferMatrix
using CairoMakie

# A lossless uniaxial crystal (calcite-like), 0.5 μm thick, in air.
# euler = (α, β, 0): β tilts the optic axis off z, α swings it azimuthally.
# (π/6, π/4, 0) tilts it OUT of the x–z plane of incidence, so p and s couple.
no, ne = 1.658, 1.486
λ = 1.0
air = Layer(λ -> 1.0, 1.0)
slab(α, β) = Layer(λ -> no, λ -> no, λ -> ne, 0.5; euler = (α, β, 0.0))

# Panel (a): sweep incidence angle with the optic axis fixed out of plane.
θdeg = 0.0:1.0:85.0
α0, β0 = π/6, π/4
res_θ = [transfer(λ, [air, slab(α0, β0), air]; θ = deg2rad(θ)) for θ in θdeg]
Rpp = getfield.(res_θ, :Rpp)
Tpp = getfield.(res_θ, :Tpp)
Rps = getfield.(res_θ, :Rps)
Tps = getfield.(res_θ, :Tps)
naive    = Rpp .+ Tpp                    # forgets both converted channels
complete = Rpp .+ Rps .+ Tpp .+ Tps      # = 1

println("Panel (a): rotated slab euler=(π/6, π/4, 0)")
println("  naive  Rpp+Tpp:      min=", round(minimum(naive), digits=5),
        "  max=", round(maximum(naive), digits=5), "  (should dip below 1)")
println("  complete Rpp+Rps+Tpp+Tps: max |budget-1| = ", maximum(abs.(complete .- 1)))

# Panel (b): sweep the optic-axis azimuth at fixed incidence; watch conversion switch on.
αdeg = 0.0:1.0:90.0
θ_fixed = deg2rad(30)
res_α = [transfer(λ, [air, slab(deg2rad(α), β0), air]; θ = θ_fixed) for α in αdeg]
Rps_α = getfield.(res_α, :Rps)
Tps_α = getfield.(res_α, :Tps)
conv_α = Rps_α .+ Tps_α
println("Panel (b): azimuth sweep at θ=30°")
println("  Rps+Tps(α=0°, optic axis in plane) = ", round(conv_α[1], digits=8), "  (≈ 0: no conversion)")
println("  Rps+Tps peaks at α≈", round(αdeg[argmax(conv_α)]), "°  =", round(maximum(conv_α), digits=5))

# --- Figure ---------------------------------------------------------------
fig = Figure(size = (980, 430))

axa = Axis(fig[1, 1];
    xlabel = "Incidence angle (°)", ylabel = "energy budget  R + T",
    title = "A lossless slab that seems to lose ~1%")
limits!(axa, 0, 85, 0.96, 1.005)
band!(axa, θdeg, naive, complete; color = (:firebrick, 0.18))
lines!(axa, θdeg, complete; color = :black, linewidth = 2.5,
    label = "Rpp + Rps + Tpp + Tps  (complete) = 1")
lines!(axa, θdeg, naive; color = :dodgerblue3, linewidth = 2.5,
    label = "Rpp + Tpp  (naive)")
text!(axa, 42, 0.978; text = "shaded gap = Rps + Tps\n(converted p→s power)",
    align = (:center, :center), fontsize = 11, color = :firebrick4)
axislegend(axa; position = :rb, framevisible = true)

axb = Axis(fig[1, 2];
    xlabel = "Optic-axis azimuth α (°)", ylabel = "converted power  (p input)",
    title = "Conversion needs an out-of-plane optic axis  (θ = 30°)")
lines!(axb, αdeg, conv_α; color = :firebrick3, linewidth = 2.5, label = "Rps + Tps")
lines!(axb, αdeg, Tps_α; color = :darkorange2, linewidth = 1.8, linestyle = :dash, label = "Tps")
lines!(axb, αdeg, Rps_α; color = :purple3, linewidth = 1.8, linestyle = :dot, label = "Rps")
scatter!(axb, [0], [conv_α[1]]; color = :black, markersize = 9)
text!(axb, 3, conv_α[1]; text = "axis in plane of incidence → no conversion",
    align = (:left, :bottom), fontsize = 11, color = :black)
limits!(axb, -2, 90, -maximum(conv_α) * 0.06, maximum(conv_α) * 1.18)
axislegend(axb; position = :rt, framevisible = true)

outpath = joinpath(@__DIR__, "..", "docs", "src", "assets", "examples", "anisotropic_energy_budget.png")
mkpath(dirname(outpath))
save(outpath, fig)
println("saved ", outpath)
