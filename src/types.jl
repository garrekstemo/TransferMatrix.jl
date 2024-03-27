struct Poynting
    out_p::SVector{3, Float64}
    in_p::SVector{3, Float64}
    out_s::SVector{3, Float64}
    in_s::SVector{3, Float64}
    refl_p::SVector{3, Float64}
    refl_s::SVector{3, Float64}

    function Poynting(out_p::T, in_p::T, out_s::T, in_s::T, refl_p::T, refl_s::T) where {T<:SVector{3, Float64}}
        new(out_p, in_p, out_s, in_s, refl_p, refl_s)
    end
end

# struct AngleResolvedResult
#     Rpp::Matrix{Float64}
#     Rss::Matrix{Float64}
#     Tpp::Matrix{Float64}
#     Tss::Matrix{Float64}
#     Γ::Vector{Matrix{ComplexF64}}
#     ξ::Matrix{ComplexF64}
# end


struct ElectricField
    z::Array{Float64}
    p::Matrix{ComplexF64}
    s::Matrix{ComplexF64}
    boundaries::Array{Float64}
end