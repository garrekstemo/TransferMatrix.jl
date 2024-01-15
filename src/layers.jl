"""
    dielectric_constant(n::Real, κ::Real)

Return the complex dielectric function from
the real and imaginary parts of the index of refraction.

The complex index of refraction, given by

        n' = n + iκ
    
(in terms of n and κ), can be used to
obtain the frequency-dependent complex dielectric function

        ε_r(ω) = ε' + iε''

via the relation

        (n + iκ)^2 = ε' + iε''.
"""
dielectric_constant(n::Real, κ::Real) = (n + κ * im)^2

"""
    dielectric_constant(n::Complex)

Return the complex dielectric function from
the complex index of refraction.
"""
dielectric_constant(n::Complex) = n^2

"""
    dielectric_constant(layer::Layer)

Return a complex dielectric function from
the index of refraction in a `Layer`` type.
"""
dielectric_constant(layer::Layer) = @. (layer.n + layer.κ * im)^2

"""
    dielectric_tensor(ε1, ε2, ε3)

Return the diagonal complex dielectric tensor

```math
\\varepsilon = 
\\begin{pmatrix}
    \\varepsilon_1 & 0 & 0 \\\\\
    0 & \\varepsilon_2  & 0 \\\\\
    0 & 0 & \\varepsilon_3
\\end{pmatrix}
```

"""
dielectric_tensor(ε1, ε2, ε3) = Diagonal(SVector{3, ComplexF64}(ε1, ε2, ε3))

"""
    permeability_tensor(μ1, μ2, μ3)

This produces the diagonal permeability tensor, 
which is identical to the way we build the `dielectric_tensor`,
and we include this function simply for completeness.
"""
permeability_tensor(μ1, μ2, μ3) = Diagonal(SVector{3, ComplexF64}(μ1, μ2, μ3))