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
We fix the wavelength of incident light and vary the angle of incidence when setting up
our `Structure`.
(It is just as simple to fix the incidence angle and calculate the transfer matrix as a function of the field wavelength.)

We start by making a `Layer` type of air and a `Layer` of glass. We'll do this for 
a wavelength of 1 μm. Since there are only two layers and the transfer matrix method
treats the first and last layers as semi-infinite, there is no need to provide a thickness
for our glass and air layers. From the examples below, you can see that there are fields for

- the material name
- the layer thickness
- a list of wavelengths
- the real part of the refractive index (corresponding to the wavelength)
- the imaginary part of the refractive index

Details about different ways to make a layer are further on in the tutorial.

```@setup tutorial
using Pkg
Pkg.add("Peaks")
Pkg.add("TransferMatrix")
Pkg.add("CairoMakie")
```

```@example tutorial
using TransferMatrix

air = Layer("Air", 0.0, [1.0e-6], [1.0], [0.0])
glass = Layer("Glass", 0.0, [1.0e-6], [1.5], [0.0])
```
### Assembling layers into a structure

Now that we have our glass and air layers, we need to assemble them into a structure
and provide the angles of the field with respect to the surface of the structure. We
do this with the `Structure` type.

```@example tutorial
θs = range(0.0, 85.0, length = 500)
s = Structure([air, glass], [1e-6], collect(θs) .* π/180)
```

The first argument is just a list of layers. The second argument is a list
of desired wavelengths. Often the refractive index data we have for 
two materials are not given for exactly the same wavelengths. 
TransferMatrix.jl uses an interpolation function to normalize
the wavelengths and complex refractive indices for all layers from 
this user-provided list of wavelengths. (Be careful if the range you
provide goes beyond the range of the data that you have!)
Now we can evaluate the angle-resolved spectrum using the function
`angle_resolved()`.

```@example tutorial
res = angle_resolved(s)
```
Let's also plot the result using the [Makie.jl](https://makie.juliaplots.org/) data visualization package.

```@example tutorial
using CairoMakie

brewster = atan(1.5) * 180 / π

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Incidence Angle (°)", ylabel = "Reflectance / Transmittance")

lines!(θs, res.Tss[:, 1], color = :firebrick4, label = "Ts")
lines!(θs, res.Tpp[:, 1], color = :orangered3, label = "Tp")
lines!(θs, res.Rss[:, 1], color = :dodgerblue4, label = "Rs")
lines!(θs, res.Rpp[:, 1], color = :dodgerblue1, label = "Rp")
vlines!(brewster, color = :dodgerblue1, linestyle = :dash)
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

## Defining a Layer

The `Layer` type is immutable. Once you make one, you can't change any of its 
characteristics later. There are two ways to define a `Layer`. The first is by 
directly filling in the argument fields:

```julia
Layer(material, thickness, λs, ns, κs)
```

You may choose any units you like for thickness and wavelength, 
but they must be the same (e.g. both are in nanometers).
The material is just a `String`. Wavelength, and the 
refractive index arguments must be Arrays, even if they just contain 
with a single item. A very manual way to make a `Layer` might look like this:

```julia
glass = Layer("Glass", 20.0e-6, [1.0e-6, 1.1e-6, 1.3e-6], [0.0, 0.0, 0.0])
```

This way works well for simple `Layers` or when you just need a single frequency
for an angle-resolved calculation, but this is a lot more work if the refractive index 
is pulled from a database or a file. The website [refractiveindex.info](https://refractiveindex.info) contains a large database of refractive indices from peer-reviewed papers.
TransferMatrix.jl uses [RefractiveIndex.jl](https://github.com/stillyslalom/RefractiveIndex.jl/tree/master) (a Julia interface to refractiveindex.info) to load refractive index data and return a `Layer` type.
Note that refractiveindex.info stores the wavelength in units of micrometers.

In the following example we load gold and assign it a thickness of 20 nm.
The optional `wavelength_unit` keyword argument represents the desired unit conversion.
In this case, if we want everything to be in meters we
must multiply each wavelength by a factor ``1e-6``.
Alternatively, we can change the thickness units to micrometers.
Whatever you do, as long as the units are consistent, the transfer matrix calculation
will perform correctly.
The `Layer` properties can be independently queried. To get the material name or thickness, 
for example, you can just type:

```@repl
audata = RefractiveMaterial("main", "Au", "Rakic-LD")
au = load_refractive_data(audata, 10e-3)
au.material
au.thickness
```

The same can be done to get wavelengths (`au.λ`) (get λ by typing `\lambda`), real refractive index (`au.n`), and complex refractive index (`au.κ`) (get κ by typing `\kappa`).


## A simple multi-layered structure

Now that we can make `Layer`s, we can put them together to calculate 
the transfer matrix for a multi-layered structure and plot the reflectance and transmittance spectra.
An important structure in optics is the Fabry-Pérot etalon, made with two parallel planar mirrors with a gap between them.
If we set the optical path length to be an integer multiple of half the wavelength, we get constructive interference and a resonance in the transmittance spectrum.
In the frequency domain, the cavity modes are evenly spaced by the free spectral range (FSR).
We will scan in the mid infrared between 4 and 6 μm and use data generated
via the [Lorentz-Drude model](https://en.wikipedia.org/wiki/Lorentz_oscillator_model) for each 10 nm-thick gold mirror. (Note that we stay in units of micrometers for the wavelength.)

```@example tutorial

λs = range(4.0, 6.0, length = 1000)
frequencies = 10^4 ./ λs  # in cm⁻¹
air = Layer("Air", 10, collect(λs), fill(1.0, length(λs)), zeros(length(λs)))
audata = RefractiveMaterial("main", "Au", "Rakic-LD")
au = load_refractive_data(audata, 10e-3)

s = Structure([air, au, air, au, air], collect(λs), [0.0])
Tp, Ts, Rp, Rs = calculate_tr(s)

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
λ_center = 5e-6
λs = collect(range(4.8, 5.2, length = 300)) .* 1e-6
θs = collect(range(0, 30, length = 100))  # degrees
νs = 10^-2 ./ λs
n1 = 2.1  
n2 = 1.6
n_bg = 1.4  # Background refractive index

# Layer thicknesses
t1 = λ_center / (4 * n1[1])
t2 = λ_center / (4 * n2[1])
t_cav = 1 * λ_center  / n_bg + 0.1e-6  # Slightly offset the cavity length to get negative detuning
```

Now we can define the absorbing material.
The amplitude of the resonance is arbitrary and is set to 5000.
The width of the resonance is set to 5 cm^-1.
The real and imaginary parts of the dielectric function are used to calculate the refractive index `n` and extinction coefficient `k`, which are used to define the `lorentzian` layer.
Since the dielectric function is unitless, we can perform these calculations in reciprocal centimeters without affecting the transfer matrix calculation, which will be done in units of micrometers.

```@example tutorial
A_0 = 5000.0
ω_0 = 10^-2 / λ_center
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
air = Layer("Air", 0.5e-6, λs, fill(1.0, length(λs)), zeros(length(λs)))
material1 = Layer("Material 1", t1, λs, fill(n1, length(λs)), zeros(length(λs)))
material2 = Layer("Material 2", t2, λs, fill(n2, length(λs)), zeros(length(λs)))
lorentzian = Layer("Lorentz oscillator", t_cav, λs, n_medium, k_medium)
nperiods = 6
unit = [material1, material2]
layers = [air, repeat(unit, nperiods)..., lorentzian, repeat(reverse(unit), nperiods)..., air]
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

s = Structure(layers, λs, deg2rad.(θs))
res = angle_resolved(s)
angle_idx, peaks = find_resonance(s, 1e-2)

θ_plot = round(rad2deg(s.θ[angle_idx]), digits=1)
T_plot = res.Tpp[angle_idx, :]
field1 = electric_field(s, s.λ[peaks[1]])
field2 = electric_field(s, s.λ[peaks[2]])

fig = Figure()
ax = Axis(fig[1, 1], title="Polariton dispersion", xlabel = "Incidence angle (°)", ylabel = "Frequency (cm⁻¹)")
heatmap!(θs, νs, res.Tpp, colormap = :deep)
fig
```

```@example tutorial

fig = Figure()

ax = Axis(fig[1, 1], title="Normal mode splitting (θ ≈ $θ_plot °)", xlabel = "Frequency  (cm⁻¹)", ylabel = "Transmittance")
lines!(νs, 4π .* k_medium .* νs ./ 1e4, color = :firebrick3, linestyle = :dash, label = "Abs")
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
lines!(field1.z .* 1e6, real(field1.p[1, :]), label = "LP")
lines!(field2.z .* 1e6, real(field2.p[1, :]), label = "UP")
hlines!(0, color = :black, linestyle = :dash, linewidth = 0.5)

ax2 = Axis(fig[2, 1], xlabel = "Distance (μm)", ylabel = "Refractive index", xticks = LinearTicks(7), yticks = [n1, n2])

# Draw DBR cavity structure
refractive_indices = [layer.n[1] for layer in layers[2:end-1]]
thicknesses = [layer.thickness for layer in layers[2:end-1]] .* 1e6

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

Here is the code to reproduce this example in its entirety:

```julia
using TransferMatrix
using Peaks
using CairoMakie

function dielectric_real(ω, p)
    A, ω_0, Γ = p
    return @. A * (ω_0^2 - ω^2) / ((ω^2 - ω_0^2)^2 + (Γ * ω)^2)
end
function dielectric_imag(ω, p)
    A, ω_0, Γ = p
    return @. A * Γ * ω / ((ω^2 - ω_0^2)^2 + (Γ * ω)^2)
end

function find_resonance(s, atol=1e-3)
    for i in eachindex(s.θ)
        pks, vals = findmaxima(res.Tpp[i, :])
        if length(pks) == 2 && isapprox(vals[1], vals[2], atol=atol)
            return i, pks
        end
    end
end

λ_center = 5e-6
λs = collect(range(4.8, 5.2, length = 300)) .* 1e-6
θs = collect(range(0, 30, length = 100))  # degrees
νs = 10^-2 ./ λs
n1 = 2.1  
n2 = 1.6
n_bg = 1.4  # Background refractive index

# Layer thicknesses
t1 = λ_center / (4 * n1[1])
t2 = λ_center / (4 * n2[1])
t_cav = 1 * λ_center  / n_bg + 0.1e-6  # Slightly offset the cavity length to get negative detuning

# Define the dielectric function inside the cavity
A_0 = 5000.0
ω_0 = 10^-2 / λ_center
Γ_0 = 5
p0 = [A_0, ω_0, Γ_0]
ε1 = dielectric_real(νs, p0) .+ n_bg^2
ε2 = dielectric_imag(νs, p0)
n_medium = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) + ε1) / 2)
k_medium = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) - ε1) / 2)

# Define the layers
air = Layer("Air", 0.5e-6, λs, fill(1.0, length(λs)), zeros(length(λs)))
material1 = Layer("Material 1", t1, λs, fill(n1, length(λs)), zeros(length(λs)))
material2 = Layer("Material 2", t2, λs, fill(n2, length(λs)), zeros(length(λs)))
lorentzian = Layer("Lorentz oscillator", t_cav, λs, n_medium, k_medium)
nperiods = 6
unit = [material1, material2]
layers = [air, repeat(unit, nperiods)..., lorentzian, repeat(reverse(unit), nperiods)..., air]

# Calculate the transmittance, reflectance
s = Structure(layers, λs, deg2rad.(θs))
res = angle_resolved(s)
angle_idx, peaks = find_resonance(s, 1e-2)

θ_plot = round(rad2deg(s.θ[angle_idx]), digits=1)
T_plot = res.Tpp[angle_idx, :]
field1 = electric_field(s, s.λ[peaks[1]])
field2 = electric_field(s, s.λ[peaks[2]])


fig = Figure(size = (400, 1000))

ax1 = Axis(fig[1, 1], title="Polariton dispersion", xlabel = "Incidence angle (°)", ylabel = "Frequency (cm⁻¹)")
heatmap!(θs, νs, res.Tpp, colormap = :deep)

ax2 = Axis(fig[2, 1], title="Normal mode splitting (θ ≈ $θ_plot °)", xlabel = "Frequency  (cm⁻¹)", ylabel = "Transmittance")
lines!(νs, 4π .* k_medium .* νs ./ 1e4, color = :firebrick3, linestyle = :dash, label = "Abs")
lines!(νs, T_plot, label = "T")
scatter!(νs[peaks], T_plot[peaks], color = :red, marker = 'x', markersize = 15)

ax3 = Axis(fig[3, 1], ylabel = "Photon field")
lines!(field1.z .* 1e6, real(field1.p[1, :]), label = "LP")
lines!(field2.z .* 1e6, real(field2.p[1, :]), label = "UP")
hlines!(0, color = :black, linestyle = :dash, linewidth = 0.5)

ax4 = Axis(fig[4, 1], xlabel = "Distance (μm)", ylabel = "Refractive index", xticks = LinearTicks(7), yticks = [n1, n2])

# Draw DBR cavity structure
refractive_indices = [layer.n[1] for layer in layers[2:end-1]]
thicknesses = [layer.thickness for layer in layers[2:end-1]] .* 1e6

prev_x = 0
prev_n = refractive_indices[1]
for (i, n) in enumerate(refractive_indices)
    current_x = sum(thicknesses[1:i])
    lines!(ax4, [prev_x, current_x], [n, n], color = :black, linewidth = 0.5)  # Plot the horizontal line

    if i > 1
        lines!(ax4, [prev_x, prev_x], [prev_n, n], color = :black, linewidth = 0.5)  # Plot the vertical line
    end
    prev_x = current_x
    prev_n = n
end

# Draw a rectangle in the cavity region
middle_i = sum(thicknesses[1:length(thicknesses) ÷ 2])
middle_f = sum(thicknesses[1:length(thicknesses) ÷ 2 + 1])
poly!(ax4, Point2f[(middle_i, n_bg), (middle_f, n_bg), (middle_f, n1), (middle_i, n1)], color = (:deepskyblue2, 0.2))

# Formatting
ylims!(ax1, 1950, 2070)
axislegend(ax2)
axislegend(ax3)
hidexdecorations!(ax3)
hideydecorations!(ax3, label = false, ticks = false, ticklabels = false)
hidexdecorations!(ax4, label = false, ticks = false, ticklabels = false)
hideydecorations!(ax4, label = false)
rowgap!(fig.layout, 3, 0)

fig
```

### Minor implementation details

Behind the scenes of the electric field calculation, a new `Structure` is actually being created and initialized for the user-supplied wavelength. The transfer matrix is again calculated for the provided wavelength and light incidence angle.
Some solution of this sort is necessary because the *exact* wavelength that you may wish to calculate for may not be part of the original `Structure`.
