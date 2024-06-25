"""
    fresnel(θ, n1, n2)

Calculate the reflectance for s-polarized and p-polarized light
given the incidence angle `θ` and indices of refraction 
of two media `n1` and `n2` at a plane interface.

[Fresnel equations](https://en.wikipedia.org/wiki/Fresnel_equations)
"""
function fresnel(θ, n1, n2)

    α = √(1 - (n1 * sin(θ) / n2)^2) / cos(θ)
    β = n2 / n1

    Rs = abs2((1/β - α) / (1/β + α))
    Rp = abs2((α - β) / (α + β))

    return Rs, Rp
end

"""
    stopband(n1, n2)

Calculate the frequency bandwidth Δf of the photonic stopband
for a distributed bragg reflector (DBR) with two alternating
materials of refractive indices `n1` and `n2`.

```math
    \\frac{\\Delta f_0}{f_0} = \\frac{4}{\\pi} \\arcsin \\left( \\frac{n_2 - n_1}{n_2 + n_1} \\right)
```

[Distributed Bragg reflector](https://en.wikipedia.org/wiki/Distributed_Bragg_reflector)
"""
function stopband(n1, n2)
    4 * asin(abs(n2 - n1) / (n2 + n1)) / π
end

"""
    dbr_reflectivity(no, ns, n1, n2, N)

Approximate the reflectivity of a DBR structure with originating medium with refractive index `no`,
substrate with index `ns`, and alternating materials with indices `n1` and `n2` and number of repetitions
`N`. The repeated pair of materials are assumed to have quarter-wave thickness ``nd = \\lambda / 4``,
where ``n`` is the refractive index, ``d`` is the layer thickness, and ``\\lambda`` is the wavelength of the light.

[Distributed Bragg reflector](https://en.wikipedia.org/wiki/Distributed_Bragg_reflector)
"""
function dbr_reflectivity(no, ns, n1, n2, N)
    r = (no * n2^(2 * N) - ns * n1^(2 * N)) / (no * n2^(2 * N) + ns * n1^(2 * N))
    return r^2
end