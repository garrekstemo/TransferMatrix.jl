module TransferMatrix

export Layer, Structure, Poynting,
       dielectric_constant,
       dielectric_tensor,
       calculate_Γ_S,
       calculate_tr,
       angle_resolved,
       electric_field,
       printstruct,
       initialize,
       read_refractive,
       load_from_yaml

using CSV
using DataInterpolations
using LinearAlgebra
import YAML

include("types.jl")
include("dataio.jl")
include("core.jl")

const ε_0 = 8.8541878128e-12
const μ_0 = 1.25663706212e-6
const c_0 = 299792458

end # module