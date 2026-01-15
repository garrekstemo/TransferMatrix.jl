# Quick Start

If you just want to get started with a transfer matrix
calculation and plot the transmittance or reflectance 
spectrum, this is the place to start.
First install TransferMatrix.jl by
typing `]` from the Julia REPL to enter package mode,
then enter

```
pkg> add TransferMatrix
```

## Quarter-wave mirror

Let's make a simple quarter-wave mirror, or
[distributed bragg reflector](https://en.wikipedia.org/wiki/Distributed_Bragg_reflector) (DBR). We will do this via alternating
layers of titanium dioxide (n = 2.13) and silica (n = 1.46) optimized
for a wavelength of 630 nm.
We'll make three periods of these two layers and and layer of air.
This 4x4 transfer matrix method simultaneously calculates
the transmittance and reflectance for s-polarized and p-polarized 
radiation.
(We are using the powerful plotting library [Makie.jl](https://makie.juliaplots.org/) to produce the figures.)


```@example dbr
using TransferMatrix
using CairoMakie
using RefractiveIndex

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
```

Note: `Layer` is parametric as `Layer{F,T}`. When annotating collections, prefer
`AbstractVector{<:Layer}` rather than `Vector{Layer}`.

Now let's try a few more periods and plot them all together
to see how the reflectance changes with increasing number of layers.
Notice that we are adding new layers directly to the structure and
not creating a new structure.

```@example dbr
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
```
