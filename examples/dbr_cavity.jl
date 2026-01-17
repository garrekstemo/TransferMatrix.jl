# DBR cavity example

using RefractiveIndex
using TransferMatrix
using GLMakie

n_air = RefractiveMaterial("other", "air", "Ciddor")
n_tio2 = RefractiveMaterial("main", "TiO2", "Sarkar")
n_sio2 = RefractiveMaterial("main", "SiO2", "Rodriguez-de_Marcos")

λ_0 = 1.0   # μm
t_tio2 = λ_0 / (4 * n_tio2(λ_0))
t_sio2 = λ_0 / (4 * n_sio2(λ_0))
t_middle = λ_0 / 2

air = Layer(n_air, t_middle)
tio2 = Layer(n_tio2, t_tio2)
sio2 = Layer(n_sio2, t_sio2)
dbr_unit = [tio2, sio2]
nperiods = 6
layers = [air, repeat(dbr_unit, nperiods)..., air, repeat(reverse(dbr_unit), nperiods)..., air];

λs = 0.85:0.001:1.2
νs = 10^4 ./ λs
Tpp = Float64[]
Tss = Float64[]
Rpp = Float64[]
Rss = Float64[]
for λ in λs
    Tpp_, Tss_, Rpp_, Rss_ = transfer(λ, layers)
    push!(Tpp, Tpp_)
    push!(Tss, Tss_)
    push!(Rpp, Rpp_)
    push!(Rss, Rss_)
end


λ = 1.0
field = efield(λ_0, layers)
ns = [real(layer.dispersion(λ_0)) for layer in layers]
bounds = [-layers[1].thickness, field.boundaries..., field.boundaries[end] + layers[end].thickness]
pushfirst!(ns, ns[1])  # Add the refractive index of the first layer for negative z


fig = Figure(size = (450, 800))
DataInspector()

Label(fig[0, :], "ZnS / MgF₂ quarter-wave stack with 3 layers", fontsize = 18)

ax1 = Axis(fig[1, 1], xlabel = "Wavenumber (cm⁻¹)", ylabel = "Transmittance / Reflectance",
            yticks = LinearTicks(5),
            xticks = LinearTicks(5),
)
lines!(νs, Tpp, label = "T")
lines!(νs, Rpp, label = "R")


# Electric field profile
ax2 = Axis(fig[2, 1], xlabel = "z position (nm)", ylabel = "Electric field")

lines!(field.z, real(field.p[1, :]))
text!(-1, 5.2, text="λ = $(Int(λ)) μm")

# Refractive index profile
ax3 = Axis(fig[3, 1], xlabel = "z position (μm)", ylabel = "Refractive index",
    xticks = 0:5,
    yticks = round.(unique(ns), digits=2),
)
stairs!(bounds, ns, color = :gray40)

axislegend(ax1)
linkxaxes!(ax2, ax3)
rowgap!(fig.layout, 3, 0.0)
hidexdecorations!(ax2, grid = false)
fig
