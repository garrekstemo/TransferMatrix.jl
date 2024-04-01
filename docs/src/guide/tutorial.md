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

We start by making a `Layer` type of air and a `Layer` of glass. We'll do this for 
a wavelength of 1 μm. Since there are only two layers and the transfer matrix method
treats the first and last layers as semi-infinite, there is no need to provide a thickness
for our glass and air layers. From the examples below, you can see that there are fields for

- the material
- the layer thickness

The material is a `RefractiveMaterial` from the RefractiveIndex.jl package.

```@setup tutorial
using Pkg
Pkg.add("Peaks")
Pkg.add("RefractiveIndex")
Pkg.add("TransferMatrix")
Pkg.add("CairoMakie")
```

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
```

Let's now plot the result using the [Makie.jl](https://makie.juliaplots.org/) data visualization package.

```@example tutorial
using CairoMakie

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

fig, ax, l = lines(frequencies, Ts)
ax.xlabel = "Frequency (cm⁻¹)"
ax.ylabel = "Transmittance"
fig
```


## Polariton dispersion in a DBR cavity

Now let's try something a little more complicated.
We will simulate an absorbing material with a resonance that coincides with a cavity mode resonance at some angle.
We will calculate the system dispersion and plot the transmittance spectrum as a function of incidence angle.
Then we will plot the electric field profile within the structure.
First we need to define the dielectric function of the absorbing material.
The real and imagary parts of the dielectric function can be defined as below:

```@example tutorial
function dielectric_real(ω, p)
    A, ω_0, Γ = p
    return @. A * (ω_0^2 - ω^2) / ((ω^2 - ω_0^2)^2 + (Γ * ω)^2)
end
function dielectric_imag(ω, p)
    A, ω_0, Γ = p
    return @. A * Γ * ω / ((ω^2 - ω_0^2)^2 + (Γ * ω)^2)
end
```
The center wavelength for the absorbing material will be 5 μm, so we set the wavelength region centered around this and define quarter-wavelength materials for the DBR in terms of this wavelength.
The refractive indices for these materials are arbitrary and, for this example, are not wavelength-dependent (although they usually are).
We will make the optical path length of the cavity region slightly greater than one wavelength to get a negative detuning between the cavity resonance and the absorbing material resonance at normal incidence (i.e. ``\Delta = |\omega_c - \omega_m| < 0``). (Here we're working in SI units.)

```@example tutorial
n_air = RefractiveMaterial("other", "air", "Ciddor")
n_tio2 = RefractiveMaterial("main", "TiO2", "Kischkat")
n_sio2 = RefractiveMaterial("main", "SiO2", "Kischkat")

λ_0 = 5.0
λs = range(4.8, 5.2, length = 50)
νs = 10^4 ./ λs
θs = range(0, 30, length = 10)
```

Now we can define the absorbing material.
The amplitude of the resonance is arbitrary and is set to 5000.
The width of the resonance is set to 5 cm^-1.
The real and imaginary parts of the dielectric function are used to calculate the refractive index `n` and extinction coefficient `k`, which are used to define the `lorentzian` layer.
Since the dielectric function is unitless, we can perform these calculations in reciprocal centimeters without affecting the transfer matrix calculation, which will be done in units of micrometers.

```@example tutorial
# absorbing material
n_bg = 1.4
A_0 = 2000.0
ω_0 = 10^4 / λ_0
Γ_0 = 5
p0 = [A_0, ω_0, Γ_0]
ε1 = dielectric_real(νs, p0) .+ n_bg^2
ε2 = dielectric_imag(νs, p0)
n_medium = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) + ε1) / 2)
k_medium = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) - ε1) / 2)
```
Now we can define the layers of the system.
For the DBR mirrors, we set the number of periods of alternating materials to 6 for high reflectivity in the region of interest.
The thickness of the two outer layers (air) is arbitrary and does not affect the result.

```@example tutorial
# dbr layers
t_tio2 = λ_0 / (4 * n_tio2(λ_0))
t_sio2 = λ_0 / (4 * n_sio2(λ_0))
t_cav = 1 * λ_0  / n_bg + 0.1  # Slightly offset the cavity length to get negative detuning

air = Layer(n_air, 2.0);
tio2 = Layer(n_tio2, t_tio2);
sio2 = Layer(n_sio2, t_sio2);
absorber = Layer(λs, n_medium, k_medium, t_cav);

nperiods = 6
unit = [tio2, sio2]
layers = [air, repeat(unit, nperiods)..., absorber, repeat(reverse(unit), nperiods)..., air];

```
After we calculate the angle-dependent transmittance and reflectance, there will be a normal mode splitting as a resulting of the coupling between the cavity mode and the absorbing material resonance.
We want to find the angle at which the amplitudes of the two modes are equal, so define the function below using the `Peaks.jl` package to do this.

```@example tutorial
using Peaks

function find_resonance(s, atol=1e-3)
    for i in eachindex(s.θ)
        pks, vals = findmaxima(res.Tpp[i, :])
        if length(pks) == 2 && isapprox(vals[1], vals[2], atol=atol)
            return i, pks
        end
    end
end
```

Now we can calculate the transmittance and reflectance spectra, find the resonance angle, and calculate the electric field at both peaks.
Then we plot the dispersion curve, the transmittance spectrum at zero detuning, and the electric field profiles at the two peaks.
Finally, we plot the refractive index profile of the structure as a function of position and display this below the field profile to see the behavior of the field at different points in the structure.


```@example tutorial
res = angle_resolved(λs, deg2rad.(θs), layers)
angle_idx, peaks = find_resonance(s, 1e-2)

θ_plot = round(rad2deg(s.θ[angle_idx]), digits=1)
T_plot = res.Tpp[angle_idx, :]
field1 = electric_field(s, λs[peaks[1]])
field2 = electric_field(s, λs[peaks[2]])

fig = Figure()
ax = Axis(fig[1, 1], title="Polariton dispersion", xlabel = "Incidence angle (°)", ylabel = "Frequency (cm⁻¹)")
heatmap!(θs, νs, res.Tpp, colormap = :deep)
fig
```

```@example tutorial

fig = Figure()

ax = Axis(fig[1, 1], title="Normal mode splitting (θ ≈ $θ_plot °)", xlabel = "Frequency  (cm⁻¹)", ylabel = "Transmittance")
lines!(νs, 4π .* k_medium .* νs, color = :firebrick3, linestyle = :dash, label = "Abs")
lines!(νs, T_plot, label = "T")
scatter!(νs[peaks], T_plot[peaks], color = :red, marker = 'x', markersize = 15)

axislegend(ax)
fig
```

We will use the following function to draw the refractive index profile of the structure.

```@example tutorial

function draw_index_profile(ax, indices, thicknesses)
    prev_x = 0
    prev_n = indices[1]
    for (i, n) in enumerate(indices)
        current_x = sum(thicknesses[1:i])
        lines!(ax, [prev_x, current_x], [n, n], color = :black, linewidth = 0.5)  # Plot the horizontal line
        if i > 1
            lines!(ax, [prev_x, prev_x], [prev_n, n], color = :black, linewidth = 0.5)  # Plot the vertical line
        end
        prev_x = current_x
        prev_n = n
    end
end
```

```@example tutorial

fig = Figure()

ax1 = Axis(fig[1, 1], ylabel = "Photon field")
lines!(field1.z .* 1e3, real(field1.p[1, :]), label = "LP")
lines!(field2.z .* 1e3, real(field2.p[1, :]), label = "UP")
hlines!(0, color = :black, linestyle = :dash, linewidth = 0.5)

ax2 = Axis(fig[2, 1], xlabel = "Distance (μm)", ylabel = "Refractive index", xticks = LinearTicks(7), yticks = [n1, n2])

# Draw DBR cavity structure
refractive_indices = [layer.n[1] for layer in layers[2:end-1]]
thicknesses = [layer.thickness for layer in layers[2:end-1]] .* 1e3

draw_index_profile(ax2, refractive_indices, thicknesses)

axislegend(ax1)
hidexdecorations!(ax1)
hideydecorations!(ax1, label = false, ticks = false, ticklabels = false)
hidexdecorations!(ax2, label = false, ticks = false, ticklabels = false)
hideydecorations!(ax2, label = false)
rowgap!(fig.layout, 1, 0)

fig
```

The electric field amplitude at both peaks is different because of the oblique incidence angle.