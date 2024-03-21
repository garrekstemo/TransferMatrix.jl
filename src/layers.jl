"""
    Layer(t, n_r, n_i)

Construct a single layer with keywords:

* `t`: thickness of the layer
* `n_r`: real part of the index of refraction
* `n_i`: imaginary part of the index of refraction
"""
struct Layer{T<:Real}
    t::T
    n_r::T
    n_i::T
end

function Layer(; t, n=1.0, n_i=0.0)
    t ≥ 0 || throw(DomainError("Layer thickness must be non-negative"))
    t, n_r, n_i = promote(t, n_r, n_i)
    return Layer(t, n_r, n_i)
end

mutable struct Structure
    layers::Vector{Layer}
    λ::Vector{Float64}
    θ::Vector{Float64}

    function Structure()
        new([], [], [])
    end

    function Structure(layers::Vector{Layer}, λs::Vector{Float64}, θs::Vector{Float64})

        new_layers = Layer[]
        for layer in layers
            new_layer = TransferMatrix.interp_data(layer, λs)
            push!(new_layers, new_layer)
        end

        new(new_layers, λs, θs)
    end

    function Structure(layers::Vector{Layer}, λs::Vector{Float64})

        new_layers = Layer[]
        for layer in layers
            new_layer = TransferMatrix.interp_data(layer, λs)
            push!(new_layers, new_layer)
        end

        new(new_layers, λs, [0.0])
    end
end


"""
    dielectric_constant(n::Real, n_i::Real)

Return the complex dielectric function from
the real and imaginary parts of the index of refraction.

The complex index of refraction, given by

        n' = n + in_i
    
(in terms of n and n_i), can be used to
obtain the frequency-dependent complex dielectric function

        ε_r(ω) = ε' + iε''

via the relation

        (n + in_i)^2 = ε' + iε''.
"""
dielectric_constant(n::Real, n_i::Real) = (n + n_i * im)^2

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
dielectric_constant(layer::Layer) = @. (layer.n + layer.n_i * im)^2

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