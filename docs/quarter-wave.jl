# Quarter-wave stack example given in Yeh pg. 110
# and Table 5.1

using TransferMatrix
using GLMakie
using CairoMakie
##

file = abspath("default_config/quarter-wave.yaml")
s = load_from_yaml(file, 1e-6)
Tpp, Tss, Rpp, Rss = calculate_tr(s)

λ_field = 1e-6
field = electric_field(s, λ_field)

CairoMakie.activate!()
f = Figure(resolution = (500, 800))
# display(f)
# DataInspector(f)

ax = Axis(f[1, 1], title = "ZnS / MgF₂ quarter-wave stack with 3 layers", xlabel = "Wavelength (nm)", ylabel = "Transmittance / Reflectance",
            yticks = LinearTicks(5), 
            xticks = LinearTicks(10))

lines!(s.λ .* 1e9, Tpp, label = "T")
lines!(s.λ .* 1e9, Rpp, label = "R")
axislegend(ax, position = :rc)

ax2 = Axis(f[2, 1], title = "Electric Field at λ = $(Int(λ_field * 1e9)) nm", xlabel = "z position (nm)", ylabel = "Field (a.u.)")

lines!(ax2, field.z .* 1e9, real(field.p[1, :]))

vlines!(ax2, field.boundaries[1] * 1e9, color = :black)
vlines!(ax2, field.boundaries[2] * 1e9, color = :firebrick4, linestyle = :dash)
vlines!(ax2, field.boundaries[3] * 1e9, color = :steelblue4, linestyle = :dash)
vlines!(ax2, field.boundaries[4] * 1e9, color = :firebrick4, linestyle = :dash)
vlines!(ax2, field.boundaries[5] * 1e9, color = :steelblue4, linestyle = :dash)
vlines!(ax2, field.boundaries[6] * 1e9, color = :firebrick4, linestyle = :dash)
vlines!(ax2, field.boundaries[end] * 1e9, color = :black)


ts = 20
text!(ax2, "ZnS", position = (110, 0.12), rotation = π/2, color = :firebrick4, textsize = ts)
text!(ax2, "MgF₂", position = (250, 0.12), rotation = π/2, color = :steelblue4, textsize = ts)
text!(ax2, "ZnS", position = (400, 0.12), rotation = π/2, color = :firebrick4, textsize = ts)
text!(ax2, "MgF₂", position = (540, 0.12), rotation = π/2, color = :steelblue4, textsize = ts)
text!(ax2, "ZnS", position = (690, 0.12), rotation = π/2, color = :firebrick4, textsize = ts)
text!(ax2, "MgF₂", position = (830, 0.12), rotation = π/2, color = :steelblue4, textsize = ts)
text!(ax2, "ZnS", position = (0.87, 0.7), space = :relative, textsize = 25)
text!(ax2, "Air", position = (0.1, 0.7), space = :relative, textsize = 25)
f