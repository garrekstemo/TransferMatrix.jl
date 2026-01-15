using RefractiveIndex
using TransferMatrix
using GLMakie
 
n_air = RefractiveMaterial("other", "air", "Ciddor")
n_glass = RefractiveMaterial("glass", "soda-lime", "Rubin-clear")[1]
air = Layer(n_air, 0.1)
glass = Layer(n_glass, 0.1)
layers = [air, glass]

λ = 1.0
θs = 0.0:1:85.0
res = sweep_angle([λ], deg2rad.(θs), layers)
brewster = rad2deg(atan(n_glass(λ)))
air.dispersion(λ)
glass.dispersion(λ)

TransferMatrix.refractive_index(n_glass)(λ)

f = Figure()
ax = Axis(f[1, 1], xlabel = "Incidence Angle (°)", ylabel = "Reflectance / Transmittance")
lines!(θs, res.Tss[:, 1], label = "Ts", color = :firebrick3)
lines!(θs, res.Tpp[:, 1], label = "Tp", color = :orangered3)
lines!(θs, res.Rss[:, 1], label = "Rs", color = :dodgerblue4)
lines!(θs, res.Rpp[:, 1], label = "Rp", color = :dodgerblue1)
vlines!(ax, brewster, color = :black, linestyle = :dash)
text!("Brewster angle\n(Rp = 0)", position = (35, 0.6))

axislegend(ax)
f