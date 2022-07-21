module TransferMatrix

export Layer, Structure, Poynting,
       dielectric_constant,
       dielectric_tensor,
       calculate_Γ_S,
       tr_from_Γ,
       tr_from_poynting,
       angle_resolved,
       printstruct,
       initialize,
       read_refractive,
       add_layer!,
       delete_layer!

using CSV
using Interpolations
using LinearAlgebra

include("types.jl")
include("dataio.jl")
include("functions.jl")
include("building.jl")

end # module