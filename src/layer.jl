"""
    Layer(material, thickness)

Construct a single layer with keywords:

* `material`: refractive material containing dispersion and extinction data (if available)
* `thickness`: thickness of the layer
"""
struct Layer
    material::RefractiveMaterial
    thickness::Real
end

function Layer(; m, t)
    t ≥ 0 || throw(DomainError("Layer thickness must be non-negative"))
    return Layer(m, t)
end

"""
    retrieve_refractive_index(material, λ)

Retrieve the refractive index for a given material at a given wavelength.
"""
function retrieve_refractive_index(material::RefractiveMaterial, λ)
    ε_real = dispersion(material, λ)
    ε_imag = 0.0
    
    try
        ε_imag = extinction(material, λ)
    catch e
        if isa(e, ArgumentError)
            # Handle the specific case of ArgumentError, which indicates missing extinction data
            ε_imag = 0.0
        else
            # Re-throw the error if it's not an ArgumentError
            rethrow(e)
        end
    end
    
    return ε_real + ε_imag * im
end

"""
    dielectric_constant(n_re::Real, n_im::Real)

Return the complex dielectric function from
the real and imaginary parts of the index of refraction.

The complex index of refraction, given by

        n' = n_re + i * n_im
    
(in terms of n_re and n_im), can be used to
obtain the frequency-dependent complex dielectric function

        ε_r(ω) = ε' + iε''

via the relation

        (n_re + i * n_im)^2 = ε' + iε''.
"""
dielectric_constant(n_re::Real, n_im::Real) = (n_re + n_im * im)^2
dielectric_constant(n::Complex) = n^2

"""
    dielectric_tensor(ε1, ε2, ε3)

Return the diagonal complex dielectric tensor

```math
\\varepsilon = 
\\begin{pmatrix}
    \\varepsilon_1 & 0 & 0 \\\\\
    0 & \\varepsilon_2  & 0 \\\\\
    0 & 0 & \\varepsilon_3
\\end{pmatrix}
```
"""
dielectric_tensor(ε1, ε2, ε3) = Diagonal(SVector{3, ComplexF64}(ε1, ε2, ε3))

"""
    permeability_tensor(μ1, μ2, μ3)

This produces the diagonal permeability tensor, 
which is identical to the way we build the `dielectric_tensor`,
and we include this function simply for completeness.
"""
permeability_tensor(μ1, μ2, μ3) = Diagonal(SVector{3, ComplexF64}(μ1, μ2, μ3))

"""
    layer_matrices(ω, ξ, layer, μ)

Calculate all parameters for a single layer, particularly
the propagation matrix and dynamical matrix so that
the overall transfer matrix can be calculated.
"""
function layer_matrices(layer, λ, ξ, μ)

    ω = 2π * c_0 / λ
    n_i = retrieve_refractive_index(layer.material, λ)
    ε_i = dielectric_constant(n_i)
    ε = dielectric_tensor(ε_i, ε_i, ε_i)

    # M = construct_M(ε, μ)
    M = SMatrix{6, 6, Complex}([
            ε[1,1] ε[1,2] ε[1,3] 0 0 0;
            ε[2,1] ε[2,2] ε[2,3] 0 0 0;
            ε[3,1] ε[3,2] ε[3,3] 0 0 0;
            0 0 0 μ 0 0;
            0 0 0 0 μ 0;
            0 0 0 0 0 μ
    ])

    a = construct_a(ξ, M)
    Δ = construct_Δ(ξ, M, a)
    q, S = calculate_q(Δ, a)
    γ = calculate_γ(ξ, q, ε, μ)
    D = dynamical_matrix(ξ, q, γ, μ)
    P = propagation_matrix(ω, q)
    return D, P, γ, q
end

"""
    printstruct(layers, unit=1e9)

Print each layer and its thickness in a somewhat 
visually useful way. Change the default unit multiplier to switch
from nanometers to micrometers. This does not affect any calculations,
only what is printed to the command line when using `printstruct`.
"""
function printstruct(layers, unit=1e9)

    unitstring = "nm"
    if unit == 1e6
        unitstring = "μm"
    end

    print("\n")
    for layer in layers

        print("-"^30, "\n")
        print("    $(layer.material), d = $(layer.thickness * unit) $(unitstring)\n")

    end
    print("-"^30, "\n")
end

"""
    find_layerbounds(layers)

Find the unitful z coordinate for all layer-layer interfaces in the structure,
with the first interface starting at z = 0.
(negative z corresponds to positions inside the first layer.)
"""
function find_layerbounds(layers)

    total_thickness = 0.0
    interface_positions = Float64[]
    
    for layer in layers
        push!(interface_positions, total_thickness + layer.thickness)
        total_thickness += layer.thickness
    end

    return interface_positions, total_thickness
end

