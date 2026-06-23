# Fabry-Pérot cavity example

using RefractiveIndex
using TransferMatrix
using CairoMakie

caf2 = RefractiveMaterial("main", "CaF2", "Malitson")
au = RefractiveMaterial("main", "Au", "Rakic-LD")
air = RefractiveMaterial("other", "air", "Ciddor")

λ_0 = 5.0
t_middle = λ_0 / 2
air = Layer(air, t_middle)
au = Layer(au, 0.01)
layers = [air, au, air, au, air]

λs = range(2.0, 6.0, length = 500)
νs = 10^4 ./ λs

Tpp = Float64[]
Tss = Float64[]
Rpp = Float64[]
Rss = Float64[]
for λ in λs
    res = transfer(λ, layers)
    push!(Tpp, res.Tpp)
    push!(Tss, res.Tss)
    push!(Rpp, res.Rpp)
    push!(Rss, res.Rss)
end

λ_min = findfirst(λs .> 4.9)
λ_max = findfirst(λs .> 5.5)
peak = findmax(Tpp[λ_min:λ_max])[2] + λ_min - 1
λ = λs[peak]

field = efield(λ, layers)


fig = Figure(size = (450, 600))

ax1 = Axis(fig[1, 1],
    title = "Fabry-Pérot Cavity",
    xlabel = "Wavelength (μm)",
    ylabel = "Transmittance / Reflectance",
    yticks = LinearTicks(5),
    xticks = LinearTicks(10),
)
lines!(νs, Tpp, label = "T")
lines!(νs, Rpp, label = "R")
scatter!(
    10^4 / λs[peak],
    Tpp[peak],
    color = :red,
    markersize = 10
)

ax2 = Axis(fig[2, 1],
    title = "Electric Field at λ = $(round(Int, λ * 1e3)) nm",
    xlabel = "z position (nm)",
    ylabel = "Field intensity (a.u.)",
)

lines!(field.z .* 1e3, real(field.p[1, :]))
vlines!(field.boundaries[1], color = :black, linestyle = :dash)
vlines!(field.boundaries[end] .* 1e3, color = :black, linestyle = :dash)

axislegend(ax1, position = :rc)

outpath = joinpath(@__DIR__, "..", "docs", "src", "assets", "examples", "fabry-perot.png")
mkpath(dirname(outpath))
save(outpath, fig)
println("saved ", outpath)
