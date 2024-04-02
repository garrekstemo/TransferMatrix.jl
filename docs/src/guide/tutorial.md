# Tutorial

## Installation

TransferMatrix.jl is a part of Julia's general registry and the source code can be found at <https://github.com/garrekstemo/TransferMatrix.jl>.
From the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), enter the package manager mode mode by typing `]`.
Then just enter the following to install the package:

```
pkg> add TransferMatrix
```

## A simple calculation

To get you up and running, let's build a simple two-layer structure of air and glass
and calculate the reflectance and transmittance to visualize the [Brewster angle](https://en.wikipedia.org/wiki/Brewster%27s_angle) for p-polarized light.
We fix the wavelength of incident light and vary the angle of incidence.

We start by making a `Layer` of air and a `Layer` of glass. We'll do this for 
a wavelength of 1 μm. Since there are only two layers and the transfer matrix method
treats the first and last layers as semi-infinite, there is no need to provide a thickness
for our glass and air layers. From the examples below, you can see that there are fields for

- the material
- the layer thickness

The material is a `RefractiveMaterial` from the RefractiveIndex.jl package. You can also initize a `Layer` with data for the refractive index and extinction coefficient.


```@example tutorial
using TransferMatrix
using RefractiveIndex

n_air = RefractiveMaterial("other", "air", "Ciddor")
n_glass = RefractiveMaterial("glass", "BK7", "SCHOTT")[1]
air = Layer(n_air, 0.1)
glass = Layer(n_glass, 0.1)
layers = [air, glass]
```

Now that we have our glass and air layers, we can iterate over the angles of incidence and compute the reflectance and transmittance for light of wavelength 1 μm.

```@example tutorial
λ = 1.0
θs = 0.0:1:85.0
result = angle_resolved([λ], deg2rad.(θs), layers)
```

Let's now plot the result using the [Makie.jl](https://makie.juliaplots.org/) data visualization package.

```@example tutorial
using CairoMakie

brewster = atan(n_glass(λ)) * 180 / π

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Incidence Angle (°)", ylabel = "Reflectance / Transmittance")
lines!(θs, result.Tss[:, 1], label = "Ts", color = :firebrick3)
lines!(θs, result.Tpp[:, 1], label = "Tp", color = :orangered3)
lines!(θs, result.Rss[:, 1], label = "Rs", color = :dodgerblue4)
lines!(θs, result.Rpp[:, 1], label = "Rp", color = :dodgerblue1)
vlines!(ax, brewster, color = :black, linestyle = :dash)
text!("Brewster angle\n(Rp = 0)", position = (35, 0.6))

axislegend(ax)
fig
```

We can see that the result of the angle-resolved calculation has four solutions: the s-wave and p-wave for both the reflected and transmitted waves. And we see that the Brewster angle
is ``\arctan\left( n_\text{glass} /n_\text{air} \right) \approx 56^{\circ}``, as expected.
Simultaneous calculation of s- and p-polarized incident waves is a feature of the 
general 4x4 transfer matrix method being used. The `angle_resolved` function
will also loop through all wavelengths so that you can plot
a color plot of wavelength and angle versus transmittance (or reflectance).


## A simple multi-layered structure

Now that we can make `Layer`s, we can put them together to calculate 
the transfer matrix for a multi-layered structure and plot the reflectance and transmittance spectra.
An important structure in optics is the Fabry-Pérot etalon, made with two parallel planar mirrors with a gap between them.
If we set the optical path length to be an integer multiple of half the wavelength, we get constructive interference and a resonance in the transmittance spectrum.
In the frequency domain, the cavity modes are evenly spaced by the free spectral range (FSR).
We will scan in the mid infrared between 4 and 6 μm and use data generated
via the [Lorentz-Drude model](https://en.wikipedia.org/wiki/Lorentz_oscillator_model) for each 10 nm-thick gold mirror. (Note that we stay in units of micrometers for the wavelength.)

```@example tutorial
caf2 = RefractiveMaterial("main", "CaF2", "Malitson")
au = RefractiveMaterial("main", "Au", "Rakic-LD")
air = RefractiveMaterial("other", "air", "Ciddor")

λ_0 = 5.0
t_middle = λ_0 / 2
layers = [Layer(air, 8.0), Layer(au, 0.01), Layer(air, t_middle), Layer(au, 0.01),  Layer(air, 8.0)];

λs = range(2.0, 6.0, length = 500)
frequencies = 10^4 ./ λs

Tp = Float64[]
Ts = Float64[]
Rp = Float64[]
Rs = Float64[]
for λ in λs
    Tp_, Ts_, Rp_, Rs_ = calculate_tr(λ, layers)
    push!(Tp, Tp_)
    push!(Ts, Ts_)
    push!(Rp, Rp_)
    push!(Rs, Rs_)
end

fig, ax, l = lines(frequencies, Ts)
ax.xlabel = "Frequency (cm⁻¹)"
ax.ylabel = "Transmittance"
fig
```

More examples are available in the `examples` folder of the package source code.