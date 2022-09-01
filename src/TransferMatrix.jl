module TransferMatrix

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
       read_refractive,
       load_from_yaml

using CSV
using DataInterpolations
using LinearAlgebra
import YAML

include("types.jl")
include("core.jl")
include("dataio.jl")

const ε_0 = 8.8541878128e-12
const μ_0 = 1.25663706212e-6
const c_0 = 299792458

end # module