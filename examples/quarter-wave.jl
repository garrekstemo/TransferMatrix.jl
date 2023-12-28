# Quarter-wave stack example given in Yeh pg. 110
# and Table 5.1

using TransferMatrix
using GLMakie

file = abspath("default_config/quarter-wave.yaml")
s = load_from_yaml(file, 1e-6)
Tpp, Tss, Rpp, Rss = calculate_tr(s)

λ_field = 1e-6
field = electric_field(s, λ_field)


fig = Figure(size = (450, 600))
DataInspector()

ax1 = Axis(fig[1, 1], title = "ZnS / MgF₂ quarter-wave stack with 3 layers", xlabel = "Wavelength (nm)", ylabel = "Transmittance / Reflectance",
            yticks = LinearTicks(5), 
            xticks = LinearTicks(10))

lines!(s.λ .* 1e9, Tpp, label = "T")
lines!(s.λ .* 1e9, Rpp, label = "R")
axislegend(ax1, position = :rc)

ax2 = Axis(fig[2, 1], title = "Electric Field at λ = $(Int(λ_field * 1e9)) nm", xlabel = "z position (nm)", ylabel = "Field intensity (a.u.)")

lines!(field.z .* 1e9, real(field.p[1, :]).^2)
vlines!(field.boundaries[1], color = :black, linestyle = :dash)
vlines!(field.boundaries[end] * 1e9, color = :black, linestyle = :dash)

fig