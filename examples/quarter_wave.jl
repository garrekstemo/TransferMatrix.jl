# Quarter-wave stack example
# from documentation

using RefractiveIndex
using TransferMatrix
using CairoMakie

n_air = RefractiveMaterial("other", "air", "Ciddor")
n_tio2 = RefractiveMaterial("main", "TiO2", "Sarkar")
n_sio2 = RefractiveMaterial("main", "SiO2", "Rodriguez-de_Marcos")

λ_0 = 0.63  # μm
t_tio2 = λ_0 / (4 * n_tio2(λ_0))
t_sio2 = λ_0 / (4 * n_sio2(λ_0))

air = Layer(n_air, 0.1)
tio2 = Layer(n_tio2, t_tio2)
sio2 = Layer(n_sio2, t_sio2)

unit = [tio2, sio2]
periods = 3
layers = [air, repeat(unit, periods)...]

λs = 0.4:0.002:1.0
Rpp = [transfer(λ, layers).Rpp for λ in λs]

fig, ax, l = lines(λs .* 1e3, Rpp)
ax.xlabel = "Wavelength (nm)"
ax.ylabel = "Reflectance"

nperiods = 6

for i in 1:nperiods
    push!(layers, tio2)
    push!(layers, sio2)
    if i%3 == 0
        Rpp_i = [transfer(λ, layers).Rpp for λ in λs]
        lines!(ax, λs .* 1e3, Rpp_i, label = "$(i + 3) periods")
    end
end

axislegend(ax)
save(joinpath(@__DIR__, "quarter_wave.png"), fig)
