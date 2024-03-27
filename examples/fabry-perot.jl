using Revise
using RefractiveIndex
using TransferMatrix
using GLMakie

caf2 = RefractiveMaterial("main", "CaF2", "Malitson")
au = RefractiveMaterial("main", "Au", "Rakic-LD")
air = RefractiveMaterial("other", "air", "Ciddor")
au.name

λ_0 = 5.0
t_middle = λ_0 / 2
layers = [Layer(air, 8.0), Layer(au, 0.01), Layer(air, t_middle), Layer(au, 0.01),  Layer(air, 8.0)];

λs = range(2.0, 6.0, length = 500)
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

λ_min = findfirst(λs .> 4.9)
λ_max = findfirst(λs .> 5.5)
peak = findmax(Tpp[λ_min:λ_max])[2] + λ_min - 1

field = electric_field(λs[peak], layers)


fig = Figure(size = (450, 600))
DataInspector()

ax1 = Axis(fig[1, 1], title = "Fabry-Pérot Cavity", xlabel = "Wavelength (μm)", ylabel = "Transmittance / Reflectance",
            yticks = LinearTicks(5), 
            xticks = LinearTicks(10))

lines!(νs, Tpp, label = "T")
scatter!(10^4 / λs[peak], Tpp[peak], color = :red, markersize = 5)
# lines!(νs, Rpp, label = "R")
axislegend(ax1, position = :rc)

ax2 = Axis(fig[2, 1], title = "Electric Field at λ = $(Int(λ * 1e3)) nm", xlabel = "z position (nm)", ylabel = "Field intensity (a.u.)")

lines!(field.z .* 1e3, real(field.p[1, :]))
vlines!(field.boundaries[1], color = :black, linestyle = :dash)
vlines!(field.boundaries[end] .* 1e3, color = :black, linestyle = :dash)

fig

##
# λs = range(3.0, 6.0, length = 1000) .* 1e-6
# θs = range(0., 50, length = 500) .* π/180

# aufile = abspath("./refractive_index_data/Au_nk_0.248-6.20um_Lorentz-Drude_Rakic1998.csv")
# caf2file = abspath("./refractive_index_data/CaF2_n_0.23-9.7um_Malitson.csv")

# air = Layer("Air", 20e-6, collect(λs), fill(1.0, length(λs)), zeros(length(λs)))
# au = read_refractive(aufile, "Au", 10e-9, div=1e6)

# s = Structure([air, au, air, au, air], collect(λs), collect(θs))

# res = angle_resolved(s)
# Tpp, Tss, Rpp, Rss = calculate_tr(s)

# λ_field = 4.025e-6
# field = electric_field(s, λ_field, 0.0, numpoints = 1000)

# ##

# fig = Figure(size = (600, 1000))
# DataInspector()

# i1 = 1
# i2 = 110
# ax1 = Axis(fig[1, 1], xlabel = "Wavelength (μm)", ylabel = "Transmittance", title = "Transmittance at $(θs[i1] * 180/π)°") #and $(round(θs[i2] * 180/π, digits=1))°")

# lines!(λs .* 1e6, res.Tss[i1, :], label = "$(round(θs[i1] * 180/π, digits=1))°")
# lines!(λs .* 1e6, res.Tss[i2, :], label = "$(round(θs[i2] * 180/π, digits=1))°")
# vlines!(λ_field * 1e6, color = :orangered, linestyle = :dash)
# axislegend(ax1)


# ax2 = Axis(fig[2, 1], title = "Field at λ = $(round(λ_field * 1e6, digits=2)) μm", xlabel = "Position within cavity (μm)", ylabel = "Electric field (a. u.)",
#         xticks = LinearTicks(5))

# lines!(field.z .* 10^6, real(field.p[1, :]), label = "Eₓ p")
# vlines!(field.boundaries[1] * 1e6, color = :black, linestyle = :dash)
# vlines!(field.boundaries[end] * 1e6, color = :black, linestyle = :dash)
# axislegend(ax2)

# ax3, hm = contourf(fig[3, 1],  θs .* 180/π, 1e-2 ./ λs, res.Tss)
# cb = Colorbar(fig[3, 2], hm, label = "Transmittance")
# ax3.title = "Transmittance contour plot"
# ax3.xlabel = "Incidence angle (°)"
# ax3.ylabel = "Wavenumber (cm⁻¹)"

# fig