module TransferMatrix

using DataInterpolations
using DelimitedFiles: readdlm
using LinearAlgebra
using RefractiveIndex
using StaticArrays


export Layer, Structure, Poynting,
       angle_resolved,
       calculate_tr,
       dielectric_constant,
       dielectric_tensor,
       electric_field,
       find_layerbounds,
       initialize,
       printstruct,
       propagation_matrix,
       load_refractive_data,
       fresnel,
       stopband,
       dbr_reflectivity


export RefractiveMaterial

include("types.jl")
include("layers.jl")
include("dataio.jl")
include("general_TMM.jl")
include("optics_functions.jl")

const ε_0::Float64 = 8.8541878128e-12
const μ_0::Float64 = 1.25663706212e-6
const c_0::Float64 = 299792458

end # module