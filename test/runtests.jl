using Test
using LinearAlgebra
using RefractiveIndex
using StaticArrays
using TransferMatrix
using Aqua

const c_0 = 299792458

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(TransferMatrix)
end

include("functions.jl")
include("layer.jl")
include("types.jl")
include("optics_functions.jl")
include("integration.jl")
