mutable struct TransferMatrixResult
    
end

struct Layer
    material::String
    thickness::Float64
    λ::Vector{Float64}
    n::Vector{Float64}
    κ::Vector{Float64}

    function Layer()
        new("Air", 0.0, [600e-9], [1.0], [0.0])
    end

    function Layer(material::String, thickness::Float64, λ, n, κ)

        if any([i < 0.0 for i in λ]) || thickness < 0.0
            error("Cannot have negative wavelength or layer thickness.")

        elseif length(λ) != length(n) || length(λ) != length(κ) || length(n) != length(κ)
            error("λ, n, and κ must all be the same length.")

        else
            new(material, thickness, λ, n, κ)
        end
    end
end

mutable struct Structure
    layers::Vector{Layer}
    λ::Vector{Float64}
    θ::Vector{Float64}
    # ξ::Vector{ComplexF64}

    function Structure()
        new([], [], [])
    end

    function Structure(layers::Vector{Layer}, λs::Vector{Float64}, θs::Vector{Float64})

        new_layers = Layer[]
        for layer in layers
            new_layer = TransferMatrix.interp_data(layer, λs)
            push!(new_layers, new_layer)
        end

        # ε_0in = TransferMatrix.dielectric_constant(layers[1])
        # ξs = @. √(ε_0in) * sin(θs)

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

struct Poynting
    out_p::Vector{Float64}
    in_p::Vector{Float64}
    out_s::Vector{Float64}
    in_s::Vector{Float64}
    refl_p::Vector{Float64}
    refl_s::Vector{Float64}

    function Poynting(out_p::Vector{Float64}, in_p::Vector{Float64}, out_s::Vector{Float64}, in_s::Vector{Float64}, refl_p::Vector{Float64}, refl_s::Vector{Float64})
        new(out_p, in_p, out_s, in_s, refl_p, refl_s)
    end
end