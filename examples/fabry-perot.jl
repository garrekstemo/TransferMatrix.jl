using TransferMatrix
using GLMakie

λs = range(3.0, 6.0, length = 1000) .* 1e-6
θs = range(0., 50, length = 500) .* π/180

aufile = abspath("./refractive_index_data/Au_nk_0.248-6.20um_Lorentz-Drude_Rakic1998.csv")
caf2file = abspath("./refractive_index_data/CaF2_n_0.23-9.7um_Malitson.csv")

air = Layer("Air", 20e-6, collect(λs), fill(1.0, length(λs)), zeros(length(λs)))
au = read_refractive(aufile, "Au", 10e-9, div=1e6)

s = Structure([air, au, air, au, air], collect(λs), collect(θs))

res = angle_resolved(s)
Tpp, Tss, Rpp, Rss = calculate_tr(s)

λ_field = 4.025e-6
field = electric_field(s, λ_field, 0.0, numpoints = 1000)

##

fig = Figure(size = (600, 1000))
DataInspector()

i1 = 1
i2 = 110
ax1 = Axis(fig[1, 1], xlabel = "Wavelength (μm)", ylabel = "Transmittance", title = "Transmittance at $(θs[i1] * 180/π)°") #and $(round(θs[i2] * 180/π, digits=1))°")

lines!(λs .* 1e6, res.Tss[i1, :], label = "$(round(θs[i1] * 180/π, digits=1))°")
lines!(λs .* 1e6, res.Tss[i2, :], label = "$(round(θs[i2] * 180/π, digits=1))°")
vlines!(λ_field * 1e6, color = :orangered, linestyle = :dash)
axislegend(ax1)


ax2 = Axis(fig[2, 1], title = "Field at λ = $(round(λ_field * 1e6, digits=2)) μm", xlabel = "Position within cavity (μm)", ylabel = "Electric field (a. u.)",
        xticks = LinearTicks(5))

lines!(field.z .* 10^6, real(field.p[1, :]), label = "Eₓ p")
vlines!(field.boundaries[1] * 1e6, color = :black, linestyle = :dash)
vlines!(field.boundaries[end] * 1e6, color = :black, linestyle = :dash)
axislegend(ax2)

ax3, hm = contourf(fig[3, 1],  θs .* 180/π, 1e-2 ./ λs, res.Tss)
cb = Colorbar(fig[3, 2], hm, label = "Transmittance")
ax3.title = "Transmittance contour plot"
ax3.xlabel = "Incidence angle (°)"
ax3.ylabel = "Wavenumber (cm⁻¹)"

fig