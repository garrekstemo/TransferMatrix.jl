# Fabry-Pérot cavity example

using Revise
using RefractiveIndex
using TransferMatrix
using GLMakie

caf2 = RefractiveMaterial("main", "CaF2", "Malitson")
au = RefractiveMaterial("main", "Au", "Rakic-LD")
air = RefractiveMaterial("other", "air", "Ciddor")

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
λ = λs[peak]

field = electric_field(λ, layers)


fig = Figure(size = (450, 600))
DataInspector()

ax1 = Axis(fig[1, 1], title = "Fabry-Pérot Cavity", xlabel = "Wavelength (μm)", ylabel = "Transmittance / Reflectance",
            yticks = LinearTicks(5), 
            xticks = LinearTicks(10))

lines!(νs, Tpp, label = "T")
scatter!(10^4 / λs[peak], Tpp[peak], color = :red, markersize = 5)
# lines!(νs, Rpp, label = "R")
axislegend(ax1, position = :rc)

ax2 = Axis(fig[2, 1], title = "Electric Field at λ = $(round(Int, λ * 1e3)) nm", xlabel = "z position (nm)", ylabel = "Field intensity (a.u.)")

lines!(field.z .* 1e3, real(field.p[1, :]))
vlines!(field.boundaries[1], color = :black, linestyle = :dash)
vlines!(field.boundaries[end] .* 1e3, color = :black, linestyle = :dash)

fig