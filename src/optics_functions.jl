"""
    fresnel(θ, n1, n2)

Calculate the reflectance for s-polarized and p-polarized light
given the incidence angle `θ` (in radians) and indices of refraction
of two media `n1` and `n2` at a plane interface.

Returns `(Rs, Rp)` where `Rs` is s-polarized reflectance and `Rp` is p-polarized reflectance.

The Fresnel equations for reflectance are:

```math
R_s = \\left| \\frac{n_1 \\cos\\theta_i - n_2 \\cos\\theta_t}{n_1 \\cos\\theta_i + n_2 \\cos\\theta_t} \\right|^2
```

```math
R_p = \\left| \\frac{n_2 \\cos\\theta_i - n_1 \\cos\\theta_t}{n_2 \\cos\\theta_i + n_1 \\cos\\theta_t} \\right|^2
```

where ``\\theta_t`` is the transmitted angle given by Snell's law:
``n_1 \\sin\\theta_i = n_2 \\sin\\theta_t``.

Special cases:
- At grazing incidence (θ → π/2), returns `(1.0, 1.0)`
- For total internal reflection (when `n1 > n2` and `θ > θ_critical`), returns `(1.0, 1.0)`

[Fresnel equations](https://en.wikipedia.org/wiki/Fresnel_equations)
"""
function fresnel(θ, n1, n2)
    cosθ = cos(θ)
    sinθ = sin(θ)

    # Handle grazing incidence (cos(θ) ≈ 0)
    if abs(cosθ) < eps(Float64)
        return (1.0, 1.0)
    end

    # Check for total internal reflection
    sin_ratio = n1 * sinθ / n2
    sin_ratio_sq = sin_ratio^2

    if sin_ratio_sq > 1.0
        # Total internal reflection: all light is reflected
        return (1.0, 1.0)
    end

    cosθt = √(1 - sin_ratio_sq)  # cos of transmitted angle

    # Standard Fresnel equations
    # Rs = |( n1 cosθ - n2 cosθt ) / ( n1 cosθ + n2 cosθt )|²
    # Rp = |( n2 cosθ - n1 cosθt ) / ( n2 cosθ + n1 cosθt )|²
    n1_cosθ = n1 * cosθ
    n2_cosθt = n2 * cosθt
    n2_cosθ = n2 * cosθ
    n1_cosθt = n1 * cosθt

    rs = (n1_cosθ - n2_cosθt) / (n1_cosθ + n2_cosθt)
    rp = (n2_cosθ - n1_cosθt) / (n2_cosθ + n1_cosθt)

    Rs = abs2(rs)
    Rp = abs2(rp)

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