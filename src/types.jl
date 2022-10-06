"""
    Layer(material::String, thickness::Float64, λ::Vector{Float64}, n::Vector{Float64}, κ::Vector{Float64})

A `Layer` stores information about a single layer,
including its material name, thickness, a list of electric field wavelengths,
and the real and imaginary parts of the refractive index
associated with these wavelengths.

Initializing a `Layer` with no arguments makes a 1 μm thick
layer of Air.
"""
struct Layer
    material::String
    thickness::Float64
    λ::Vector{Float64}
    n::Vector{Float64}
    κ::Vector{Float64}

    function Layer()
        new("Air", 0.0, [1e-6], [1.0], [0.0])
    end

    function Layer(material, thickness, λ, n, κ)

        if any([i < 0.0 for i in λ]) || thickness < 0.0
            error("Cannot have negative wavelength or layer thickness.")

        elseif length(λ) != length(n) || length(λ) != length(κ) || length(n) != length(κ)
            error("λ, n, and κ must all be the same length.")

        else
            new(material, thickness, λ, n, κ)
        end
    end
end


"""
    Structure

The `Structure` is a mutable type that stores a Vector of `Layer` types, along with
a list of field wavelengths and incident angles to calculate on.
"""
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

struct TransferMatrixResult
    tm::Vector{Matrix{ComplexF64}}
    poynting::Vector{Poynting}
    ξ::Vector{Complex}
end

struct AngleResolvedResult
    Rpp::Matrix{Float64}
    Rss::Matrix{Float64}
    Tpp::Matrix{Float64}
    Tss::Matrix{Float64}
    Γ::Vector{Matrix{ComplexF64}}
    ξ::Matrix{ComplexF64}
end


struct ElectricField
    z::Array{Float64}
    p::Matrix{ComplexF64}
    s::Matrix{ComplexF64}
    boundaries::Array{Float64}
end