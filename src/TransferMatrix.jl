module TransferMatrix

using DataInterpolations
using LinearAlgebra
using PrecompileTools
using RefractiveIndex
using StaticArrays

export Layer,
       TransferResult,
       sweep_angle,
       transfer,
       dielectric_constant,
       dielectric_tensor,
       ElectricField,
       efield,
       find_bounds,
       fresnel,
       fresnel_coefficients,
       airy,
       stopband,
       dbr_reflectivity,
       refractive_index,
       sweep_thickness,
       isanisotropic,
       get_refractive_indices,
       get_euler_angles,
       isrotated,
       euler_rotation_matrix,
       rotate_dielectric_tensor

include("matrix_constructors.jl")
include("layer.jl")
include("general_TMM.jl")
include("optics_functions.jl")

const ε_0::Float64 = 8.8541878128e-12
const μ_0::Float64 = 1.25663706212e-6
const c_0::Float64 = 299792458

# Precompile common workloads to reduce time-to-first-execution
@setup_workload begin
    # Use simple constant-index dispersion functions (no RefractiveIndex.jl dependency)
    n_air = λ -> 1.0 + 0.0im
    n_film = λ -> 1.5 + 0.0im
    n_sub = λ -> 1.45 + 0.0im

    λ_0 = 1.0
    d_film = λ_0 / 6.0  # ~quarter-wave

    @compile_workload begin
        # Layer construction
        air = Layer(n_air, 0.1)
        film = Layer(n_film, d_film)
        sub = Layer(n_sub, 0.5)
        layers = [air, film, sub]

        # Core TMM calculation paths
        transfer(λ_0, layers)
        transfer(λ_0, layers; θ=0.3)

        # Electric field calculation
        efield(λ_0, layers; dz=0.01)
    end
end

end # module
