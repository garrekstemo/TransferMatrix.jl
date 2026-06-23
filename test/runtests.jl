using Test
using LinearAlgebra
using RefractiveIndex
using StaticArrays
using TransferMatrix
using Aqua
using Unitful

const c_0 = 299792458

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(TransferMatrix; deps_compat=(check_extras=false, ignore=[:LinearAlgebra],))
end

include("functions.jl")
include("layer.jl")
include("types.jl")
include("optics_functions.jl")
include("dispersion.jl")
include("integration.jl")
include("circular.jl")
include("sheets.jl")
include("refractiveindex.jl")
include("unitful.jl")
