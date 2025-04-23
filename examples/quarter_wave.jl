# Quarter-wave stack example
# from documentation

using RefractiveIndex
using TransferMatrix
using GLMakie

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
Rpp = Float64[]
for λ in λs
    Tpp_, Tss_, Rpp_, Rss_ = calculate_tr(λ, layers)
    push!(Rpp, Rpp_)
end


f, ax, l = lines(λs .* 1e3, Rpp)
ax.xlabel = "Wavelength (nm)"
ax.ylabel = "Reflectance"
f

##

nperiods = 6

for i in 1:nperiods
    push!(layers, tio2)
    push!(layers, sio2)
    Rpp = Float64[]
    if i%3 == 0
        for λ in λs
            Tpp_, Tss_, Rpp_, Rss_ = calculate_tr(λ, layers)
            push!(Rpp, Rpp_)
        end
        lines!(ax, λs .* 1e3, Rpp, label = "$(i + 3) periods")
    end
end

axislegend(ax)
f