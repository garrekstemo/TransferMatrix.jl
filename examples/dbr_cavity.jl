# DBR cavity example

using RefractiveIndex
using TransferMatrix
using GLMakie

n_air = RefractiveMaterial("other", "air", "Ciddor")
n_tio2 = RefractiveMaterial("main", "TiO2", "Sarkar")
n_sio2 = RefractiveMaterial("main", "SiO2", "Rodriguez-de_Marcos")

λ_0 = 1.0   # μm
nperiods = 6
t_tio2 = λ_0 / (4 * n_tio2(λ_0))
t_sio2 = λ_0 / (4 * n_sio2(λ_0))
t_middle = λ_0 / 2

air = Layer(n_air, t_middle)
tio2 = Layer(n_tio2, t_tio2)
sio2 = Layer(n_sio2, t_sio2)
dbr_unit = [tio2, sio2]
layers = [air, repeat(dbr_unit, nperiods)..., air, repeat(reverse(dbr_unit), nperiods)..., air];

λs = 0.8:0.001:1.2
νs = 10^4 ./ λs
Tpp = Float64[]
Tss = Float64[]
Rpp = Float64[]
Rss = Float64[]
for λ in λs
    Tpp_, Tss_, Rpp_, Rss_ = calculate_tr(λ, layers)
    push!(Tpp, Tpp_)
    push!(Tss, Tss_)
    push!(Rpp, Rpp_)
    push!(Rss, Rss_)
end


λ = 1.0
field = electric_field(λ_0, layers)


fig = Figure(size = (450, 600))
DataInspector()

ax1 = Axis(fig[1, 1], title = "ZnS / MgF₂ quarter-wave stack with 3 layers", xlabel = "Wavelength (μm)", ylabel = "Transmittance / Reflectance",
            yticks = LinearTicks(5), 
            xticks = LinearTicks(10))

lines!(λs, Tpp, label = "T")
lines!(λs, Rpp, label = "R")
axislegend(ax1, position = :rc)

ax2 = Axis(fig[2, 1], title = "Electric Field at λ = $(Int(λ * 1e3)) nm", xlabel = "z position (nm)", ylabel = "Field intensity (a.u.)")

lines!(field.z .* 1e3, real(field.p[1, :]))
vlines!(field.boundaries[1], color = :black, linestyle = :dash)
vlines!(field.boundaries[end] .* 1e3, color = :black, linestyle = :dash)

fig