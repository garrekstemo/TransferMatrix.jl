module TransferMatrix

using DataInterpolations
using LinearAlgebra
using RefractiveIndex
using StaticArrays

export Layer,
       sweep_angle,
       calculate_tr,
       dielectric_constant,
       dielectric_tensor,
       electric_field,
       find_bounds,
       fresnel,
       stopband,
       dbr_reflectivity,
       refractive_index,
       sweep_thickness

include("matrix_constructors.jl")
include("layer.jl")
include("general_TMM.jl")
include("optics_functions.jl")

const ε_0::Float64 = 8.8541878128e-12
const μ_0::Float64 = 1.25663706212e-6
const c_0::Float64 = 299792458

end # module
