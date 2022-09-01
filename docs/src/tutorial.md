# Tutorial

## Getting started

### Installation

TransferMatrix.jl is not a part of Julia's general registry, but can be found at <https://github.com/garrekstemo/TransferMatrix.jl>. Julia makes it easy to 
install unregistered packages. From the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), enter the Pkg REPL (package manager) mode by typing `]` and then use `add` to add the GitHub URL to TransferMatrix.jl.

```
pkg> add https://github.com/garrekstemo/TransferMatrix.jl
```

### A simple calculation

To get you up and running, let's build a simple two-layer structure of air and glass
and calculate the reflectance and transmittance spectrum while changing the
angle of incidence for the incoming electromagnetic field. It is just as simple to
fix the incidence angle and calculate the transfer matrix as a function of the field wavelength.

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

```julia
using TransferMatrix

air = Layer("Air", 0.0, [1.0e-6], [1.0], [0.0])
glass = Layer("Glass", 0.0, [1.0e-6], [1.5], [0.0])
```
#### Assembling layers into a structure

Now that we have our glass and air layers, we need to assemble them into a structure
and provide the angles of the field with respect to the surface of the structure. We
do this with the `Structure` type.

```julia
θs = collect(range(0., 80., length = 500)) .* π/180
s = Structure([air, glass], [1.0e-6], θs)
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

```julia
res = angle_resolved(s)
```
Let's also plot the result using the [Makie.jl](https://makie.juliaplots.org/) data visualization package.
```julia
using CairoMakie

f = Figure()
ax = Axis(f[1, 1], xlabel = "Incidence Angle (°)", ylabel = "Reflectance / Transmittance")

lines!(s.θ .* 180/π, res.Rss[:, 1], color = :dodgerblue4, label = "Rss")
lines!(s.θ .* 180/π, res.Rpp[:, 1], color = :dodgerblue1, label = "Rpp")
lines!(s.θ .* 180/π, res.Tss[:, 1], color = :firebrick4, label = "Tss")
lines!(s.θ .* 180/π, res.Tpp[:, 1], color = :orangered3, label = "Tpp")

axislegend(ax)

f
```
```@raw html
<p align="center">
    <img src="../assets/air-glass_example.svg">
</p>
```

We can see that the result of the angle-resolved calculation is Julia type
with four solutions: the s-wave and p-wave for both the reflected and transmitted waves.
Simultaneous calculation of s- and p- incident waves is a feature of the 
general 4x4 transfer matrix method being used. The `angle_resolved` function
will also loop through all wavelengths provided so that you can plot
a heatmap of wavelength and angle versus transmittance (or reflectance).

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
TransferMatrix.jl makes it easy to use a CSV file downloaded from this website with the `read_refractive()` function, which returns a `Layer` type. You can choose 
to include (or not) the real/imaginary part of the complex refractive index. Note that 
this refractiveindex.info stores wavelength information in units of micrometers.

```julia
glass = read_refractive("/path/to/file.csv", "Glass", 20e-6, div=1e6)
```

The `div` keyword argument at the end of `read_refractive` represents the 
desired unit conversion. In this case, if we want everything to be in meters we
must divide each wavelength in the raw data by ``10^6``. Alternatively, we can 
change the thickness units to micrometers. Just be careful that all `Layer`s
in your multi-layered structure have consistent units.

The `Layer` properties can be independently queried. To get the material name or thickness, 
for example, you can just type:

```julia
julia> glass.material
"Glass"
julia> glass.thickness
2.0e-5
```

The same can be done to get wavelengths (`glass.λ`) (get λ by typing `\lambda`), real refractive index (`glass.n`), and complex refractive index (`glass.κ`) (get κ by typing `\kappa`). (A very useful feature of Julia is built-in support for unicode characters).

## A simple multi-layered structure

Now that we can make `Layer`s, we can put them together to calculate 
the transfer matrix for a multi-layered structure and plot
the reflectance and transmittance spectra. A very important structure
in optics is the Fabry-Pérot etalon, which can be made with two
mirrors that face each other with a gap between them. Let's make this a 
little more challenging (and experimentally realistic) by 
mounting the mirrors onto windows that are transparent in our
chosen wavelength range. The gap will simply be 20 μm of air.
We will scan in the mid infrared between 4 and 6 μm and use data generated
via the [Lorentz-Drude model](https://en.wikipedia.org/wiki/Lorentz_oscillator_model) for each 10 nm-thick gold layer and experimental data
from Malitson, 1963 for the windows.


```julia
λs = range(4.0, 6.0, length = 1000) .* 1e-6
aufile = abspath("./refractive_index_data/Au_nk_0.248-6.20um_Lorentz-Drude_Rakic1998.csv")
caf2file = abspath("./refractive_index_data/CaF2_n_0.23-9.7um_Malitson.csv")

air = Layer("Air", 20e-6, collect(λs), fill(1.0, length(λs)), zeros(length(λs)))
au = read_refractive(aufile, "Au", 10e-9, div=1e6)
caf2 = read_refractive(caf2file, "CaF2", 5.0e-6, div=1e6)

s = Structure([caf2, au, air, au, caf2], collect(λs), [0.0])
Tpp, Tss, Rpp, Rss = calculate_tr(s)
```

```@raw html
<p align="center">
    <img src="../assets/fabry-perot_example.svg">
</p>
```

It is easy to combine manually-generated data and experimental data
to calculate the global transfer matrix.


## Periodic structures with a YAML config file

With the tools described above, it is pretty easy to make any kind of multi-layered
structure that you want. A very large structure, maybe one with repeating layers,
is also not difficult to build. However, it would be useful to be able to load a 
structure from a set of parameters so that it can be saved or even shared with others for better reproducibility.
To do this, a `Structure` can be loaded directly from a YAML file (using [YAML.jl](https://github.com/JuliaData/YAML.jl)). 

Let's demonstrate this with a [quarter-wave stack](https://en.wikipedia.org/wiki/Distributed_Bragg_reflector), which is a periodic structure
with two alternative layers where the thickness is one fourth of the wavelength
within the medium, or the relation

```math
d_i = \frac{\lambda}{4n_i},
```
where ``d_i`` is the thickness for layer ``i``, ``\lambda`` is the electric field wavelength in vacuum, and ``n_i`` is the index of refraction for layer ``i``.

Here is the YAML file that gives us this:

```
# All units in micrometers

min_wavelength: 0.5
max_wavelength: 2.0
n_wavelengths: 1000

# Theta in radians
theta_i: 0.0
theta_f: 0.0
n_theta: 1

layers:

    layer1:
        material: Air
        thickness: 0.5
        wavelength: 1.0
        refractive_index: 1.0
        extinction_coeff: 0.0

    layer2:
        periods: 3
        layers:
          layer1:
            material: "ZnS n=2.32"
            thickness: 0.1077
            refractive_filename: ZnS_n_0.4-14um_Amotchkina2020.csv
          layer2:
            material: "MgF2 n=1.36"
            thickness: 0.1838
            refractive_filename: MgF2_nk_0.03-2.0um_Rodriguez2017.csv

    layer3: 
        material: "ZnS n=2.32"
        thickness: 0.1838
        refractive_filename: ZnS_n_0.4-14um_Amotchkina2020.csv
```

At the top, we set the minimum and maximum wavelengths to be calculated, with the 
total number of points set to 1000. We won't do an angle-resolved calculation, so
initial and final ``\theta`` are zero (with `n_theta` just one). Next the multi-layered structure is defined under `layers:`. The first layer is just air,
so we don't need a file for this and can set the refractive index for a single wavelength. These values and wavelength are expanded when the YAML file is loaded.
The final layer (`layer3`) is zinc sulfide (ZnS), which is loaded from a CSV file
downloaded from [refractiveindex.info](https://refractiveindex.info). The middle section (`layer2`) is where we define the periodic structure. The first item
determines how many periods we want. In this case, `periods: 3`. Then we 
set the layers that alternate under another `layers` section. You can see that
the first layer `layer1` is zinc sulfide and the second layer `layer2` is magnesium fluoride (MgF₂). The number of times that these two layers are repeated can easily be changed by editing the `periods` field. You can even have multiple repeating structures within the global structure for arbitrarily complicated structures.

To load the YAML config file into a `Structure`, we use the `load_from_yaml()` function. You can find this example in the `default_config` folder
in the `quarter-wave.yaml` file.

```julia
s = load_from_yaml(".../default_config/quarter-wave.yaml")
Tpp, Tss, Rpp, Rss = calculate_tr(s)
```

This particular example is taken from Pochi Yeh's [*Optical Waves in Layered Media*](https://www.wiley.com/en-us/Optical+Waves+in+Layered+Media-p-9780471731924) on page 110. You can find the reflectivity of this structure in Table 5.1. The values we have calculated are slightly different since the data used here is slightly dispersive with wavelength, but the example in the book takes flat refractive index values. You can try plotting this structure for an increasing number of periods and observe how the reflectance near ``\lambda`` = 1 μm increases.

```julia
using CairoMakie

f = Figure()
ax = Axis(f[1, 1], title = "ZnS / MgF quarter-wave with 3 layers", xlabel = "Wavelength (nm)", ylabel = "Transmittance / Reflectance")

lines!(s.λ .* 1e9, Tpp, label = "T")
lines!(s.λ .* 1e9, Rpp, label = "R")
axislegend(ax, position = :rc)

f
```

```@raw html
<p align = "center">
    <img src = "/assets/quarter-wave_example.svg">
</p>
```

## Electric field

The electric field can be calculated as a function of position within the
layered structure using the `electric_field()` function, which takes
the `Structure` and desired wavelength `λ`, as well as optional argument angle of incidence `θ` and optional keyword argument for the number of data points `numpoints`. We can plot the field profile for the [distributed bragg reflector](https://www.rp-photonics.com/bragg_mirrors.html) (DBR) we constructed in the previous section. Let's do this for `λ` = 1 μm.

```julia
λ_field = 1e-6
field = electric_field(s, λ_field)

f = Figure()
ax = Axis(f[1, 1], title = "Electric Field Profile at λ = $(Int(λ_field * 1e9)) nm", xlabel = "z position (nm)", ylabel = "Field intensity (a.u.)")

lines!(field.z .* 1e9, real(field.p[1, :]).^2)

vlines!(field.boundaries[1], color = :black)
vlines!(field.boundaries[end] * 1e9, color = :black)

f
```

```@raw html
<p align = "center">
    <img src = "/assets/dbr-field_example.svg">
</p>
```

The electric field result contains the position `z` within the structure,
the (x, y, z) components (corresponding to the first, second, and third components of the Array) of the p-polarized and s-polarized light, and the positions of the
layer interfaces. So for example, to get the x-component of the p-polarized field along all of `z`, we would call `field.p[1, :]`, as we have done above.

The layer `boundaries` is useful for plotting (as shown in the figure above) and checking that the in-plane components are continuous throughout the structure, as required by [Maxwell's interface conditions](https://en.wikipedia.org/wiki/Interface_conditions_for_electromagnetic_fields).

### Minor implementation details

Behind the scenes, a new `Structure` is actually being created and initialized for the user-supplied wavelength. The transfer matrix is again calculated for the provided wavelength and light incidence angle.
Some solution of this sort is necessary because the *exact* wavelength that you may wish to calculate for may not be part of the original `Structure`.
