using TransferMatrix
using GLMakie
using CairoMakie

air = Layer("Air", 0.0, [1.0e-6], [1.0], [0.0])
glass = Layer("Glass", 0.0, [1.0e-6], [1.5], [0.0])

θs = collect(range(0., 80., length = 500)) .* π/180
s = Structure([air, glass], [1.0e-6], θs)
res = angle_resolved(s)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Incidence Angle (°)", ylabel = "Reflectance / Transmittance")

lines!(s.θ .* 180/π, res.Rss[:, 1], color = :dodgerblue4, label = "Rss")
lines!(s.θ .* 180/π, res.Rpp[:, 1], color = :dodgerblue1, label = "Rpp")
lines!(s.θ .* 180/π, res.Tss[:, 1], color = :firebrick4, label = "Tss")
lines!(s.θ .* 180/π, res.Tpp[:, 1], color = :orangered3, label = "Tpp")

axislegend(ax)

fig