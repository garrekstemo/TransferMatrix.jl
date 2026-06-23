# Energy conservation in a rotated anisotropic slab: the missing cross-polarization.
#
# A lossless dielectric slab must obey R + T = 1. But take a uniaxial crystal and
# rotate its optic axis OUT of the plane of incidence, and the naive budget
#
#     R + T  =  Rpp + Tpp  ≈  0.998
#
# falls short of 1 — even though nothing absorbs. The slab is not leaking energy:
# a rotated anisotropic crystal CONVERTS polarization. A purely p-polarized input
# comes back with a small s-polarized component (Rps), and that reflected
# cross-polarized light is exactly the piece the naive budget forgets. The
# complete, polarization-resolved budget
#
#     Rpp + Rps + Tpp  =  1
#
# closes to machine precision. (Tpp from the Poynting flux already includes the
# cross-transmitted light, so it is not added twice.)
#
# Panel (a): the naive vs. complete budget vs. incidence angle — the shaded gap IS Rps.
# Panel (b): Rps vs. the optic-axis azimuth — conversion vanishes when the axis lies
#            in the plane of incidence and grows as it tilts out.
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
naive    = Rpp .+ Tpp            # forgets the p→s reflection
complete = Rpp .+ Rps .+ Tpp     # = 1

println("Panel (a): rotated slab euler=(π/6, π/4, 0)")
println("  naive  Rpp+Tpp:      min=", round(minimum(naive), digits=5),
        "  max=", round(maximum(naive), digits=5), "  (should dip below 1)")
println("  complete Rpp+Rps+Tpp: max |budget-1| = ", maximum(abs.(complete .- 1)))

# Panel (b): sweep the optic-axis azimuth at fixed incidence; watch Rps switch on.
αdeg = 0.0:1.0:90.0
θ_fixed = deg2rad(30)
Rps_α = [transfer(λ, [air, slab(deg2rad(α), β0), air]; θ = θ_fixed).Rps for α in αdeg]
println("Panel (b): azimuth sweep at θ=30°")
println("  Rps(α=0°, optic axis in plane) = ", round(Rps_α[1], digits=8), "  (≈ 0: no conversion)")
println("  Rps peaks at α≈", round(αdeg[argmax(Rps_α)]), "°  =", round(maximum(Rps_α), digits=5))

# --- Figure ---------------------------------------------------------------
fig = Figure(size = (980, 430))

axa = Axis(fig[1, 1];
    xlabel = "Incidence angle (°)", ylabel = "energy budget  R + T",
    title = "A lossless slab that seems to lose ~0.2%")
limits!(axa, 0, 85, 0.99, 1.002)
band!(axa, θdeg, naive, complete; color = (:firebrick, 0.18))
lines!(axa, θdeg, complete; color = :black, linewidth = 2.5,
    label = "Rpp + Rps + Tpp  (complete) = 1")
lines!(axa, θdeg, naive; color = :dodgerblue3, linewidth = 2.5,
    label = "Rpp + Tpp  (naive)")
text!(axa, 42, 0.9988; text = "shaded gap = Rps\n(p→s reflection)",
    align = (:center, :center), fontsize = 11, color = :firebrick4)
axislegend(axa; position = :rb, framevisible = true)

axb = Axis(fig[1, 2];
    xlabel = "Optic-axis azimuth α (°)", ylabel = "cross-pol reflectance  Rps",
    title = "Conversion needs an out-of-plane optic axis  (θ = 30°)")
lines!(axb, αdeg, Rps_α; color = :firebrick3, linewidth = 2.5)
scatter!(axb, [0], [Rps_α[1]]; color = :black, markersize = 9)
text!(axb, 3, Rps_α[1]; text = "axis in plane of incidence → Rps = 0",
    align = (:left, :bottom), fontsize = 11, color = :black)
limits!(axb, -2, 90, -maximum(Rps_α) * 0.06, maximum(Rps_α) * 1.18)

outpath = joinpath(@__DIR__, "..", "docs", "src", "assets", "examples", "anisotropic_energy_budget.png")
mkpath(dirname(outpath))
save(outpath, fig)
println("saved ", outpath)
