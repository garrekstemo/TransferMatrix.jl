using Revise
using RefractiveIndex
using TransferMatrix
using GLMakie
 
n_air = RefractiveMaterial("other", "air", "Ciddor")
n_glass = RefractiveMaterial("glass", "BK7", "SCHOTT")[1]
n_glass(1.0)
air = Layer(n_air, 0.1)
glass = Layer(n_glass, 0.1)
layers = [air, glass]

λ = 1.0
θs = 0.0:1:85.0

Tpp = Float64[]
Tss = Float64[]
Rpp = Float64[]
Rss = Float64[]
for θ in θs
    Tpp_, Tss_, Rpp_, Rss_ = calculate_tr(λ, layers, deg2rad(θ))
    push!(Tpp, Tpp_)
    push!(Tss, Tss_)
    push!(Rpp, Rpp_)
    push!(Rss, Rss_)
end

brewster = atan(n_glass(λ)) * 180 / π

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Incidence Angle (°)", ylabel = "Reflectance / Transmittance")
lines!(θs, Tss, label = "Ts", color = :firebrick3)
lines!(θs, Tpp, label = "Tp", color = :orangered3)
lines!(θs, Rss, label = "Rs", color = :dodgerblue4)
lines!(θs, Rpp, label = "Rp", color = :dodgerblue1)
vlines!(ax, brewster, color = :black, linestyle = :dash)
text!("Brewster angle\n(Rp = 0)", position = (35, 0.6))

axislegend(ax)
fig