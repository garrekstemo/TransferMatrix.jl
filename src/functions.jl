const c_0 = 299792458

"""
Return a complex dielectric function from
the index of refraction in a Layer type.

The complex index of refraction, given by

        n* = n + iκ
    
(in terms of n and \\kappa), can be used to
obtain the frequency-dependent complex dielectric function
    
        ε_r(ω) = ε' + iε''

via the relation

        (n + iκ)^2 = ε' + iε''.
"""
function dielectric_constant(layer::Layer)
    return @. (layer.n + layer.κ*im)^2
end

"""
Return the complex dielectric function from
the real and imaginary parts of the index of refraction.
"""
function dielectric_constant(n::Float64, κ::Float64)
    return (n + κ*im)^2
end

"""
Return the complex dielectric function from
the complex index of refraction.
"""
function dielectric_constant(n::ComplexF64)
    return n^2
end

"""
Return the diagonal complex dielectric tensor

     [ε1 0  0
ε =   0  ε2 0
      0  0  ε3]
"""
function dielectric_tensor(ε1, ε2, ε3)
    [ε1 0 0 ; 0 ε2 0 ; 0 0 ε3]
end

function dielectric_tensor(layer::Layer)
    
    εs = dielectric_constant(layer)
    ε_tensor = zeros(length(εs))

    for (i, ε_i) in enumerate(εs)
        ε_tensor[i] = dielectric_tensor(ε_i, ε_i, ε_i)
    end
    ε_tensor
end


function construct_a(ξ, M)

    a = zeros(ComplexF64, 6, 6)

    b = M[3, 3] * M[6, 6] - M[3, 6] * M[6, 3]

    a[3,1] = (M[6,1] * M[3,6] - M[3,1] * M[6,6]) / b
    a[3,2] =((M[6,2] - ξ) * M[3,6] - M[3,2] * M[6,6]) / b
    a[3,4] = (M[6,4] * M[3,6] - M[3,4] * M[6,6]) / b
    a[3,5] = (M[6,5] * M[3,6] - (M[3,5] + ξ) * M[6,6]) / b

    a[6,1] = (M[6,3] * M[3,1] - M[3,3] * M[6,1]) / b
    a[6,2] = (M[6,3] * M[3,2] - M[3,3] * (M[6,2] - ξ)) / b
    a[6,4] = (M[6,3] * M[3,4] - M[3,3] * M[6,4]) / b
    a[6,5] = (M[6,3] * (M[3,5] + ξ) - M[3,3] * M[6,5]) / b
    
    return a
end

"""
Construct the reordered matrix Δ in terms of the elements of M and a
and the in-plane reduced wavevector ξ = k_x / k_0.
"""
function constructΔ(ξ, M, a)

    Δ = Array{ComplexF64}(undef, 4, 4)

    Δ[1,1] =  M[5,1] + (M[5,3] + ξ) * a[3,1] + M[5,6] * a[6,1]
    Δ[1,2] =  M[5,5] + (M[5,3] + ξ) * a[3,5] + M[5,6] * a[6,5]
    Δ[1,3] =  M[5,2] + (M[5,3] + ξ) * a[3,2] + M[5,6] * a[6,2]
    Δ[1,4] = -M[5,4] - (M[5,3] + ξ) * a[3,4] - M[5,6] * a[6,4]

    Δ[2,1] =  M[1,1] + M[1,3] * a[3,1] + M[1,6] * a[6,1]
    Δ[2,2] =  M[1,5] + M[1,3] * a[3,5] + M[1,6] * a[6,5]
    Δ[2,3] =  M[1,2] + M[1,3] * a[3,2] + M[1,6] * a[6,2]
    Δ[2,4] = -M[1,4] - M[1,3] * a[3,4] - M[1,6] * a[6,4]

    Δ[3,1] = -M[4,1] - M[4,3] * a[3,1] - M[4,6] * a[6,1]
    Δ[3,2] = -M[4,5] - M[4,3] * a[3,5] - M[4,6] * a[6,5]
    Δ[3,3] = -M[4,2] - M[4,3] * a[3,2] - M[4,6] * a[6,2]
    Δ[3,4] =  M[4,4] + M[4,3] * a[3,4] + M[4,6] * a[6,4]

    Δ[4,1] =  M[2,1] + M[2,3] * a[3,1] + (M[2,6] - ξ) * a[6,1]
    Δ[4,2] =  M[2,5] + M[2,3] * a[3,5] + (M[2,6] - ξ) * a[6,5]
    Δ[4,3] =  M[2,2] + M[2,3] * a[3,2] + (M[2,6] - ξ) * a[6,2]
    Δ[4,4] = -M[2,4] - M[2,3] * a[3,4] - (M[2,6] - ξ) * a[6,4]

    return Δ
end



"""
From Berreman, 1972, Ψ is the column matrix:

    [ Ex
Ψ =   Hy
      Ey
     -Hx ]

for a right-handed Cartesian coordinate system with
the z-axis along the normal to the multilayer structure.
"""
function poynting(Ψ, a)

    S = Array{ComplexF64}(undef, 3, 4)

    for m in 1:4
        Ex =  Ψ[1, m]
        Ey =  Ψ[3, m]
        Hx = -Ψ[4, m]
        Hy =  Ψ[2, m]

        Ez = a[3,1] * Ex + a[3,2] * Ey + a[3,4] * Hx + a[3,5] * Hy
        Hz = a[6,1] * Ex + a[6,2] * Ey + a[6,4] * Hx + a[6,5] * Hy
        
        S[1, m] = Ey * Hz - Ez * Hy
        S[2, m] = Ez * Hx - Ex * Hz
        S[3, m] = Ex * Hy - Ey * Hx
    end
    return S
end

function poynting(ξ, q_in, q_out, γ_in, γ_out, t_coefs, r_coefs)

    # create the wavevector in the first layer
    k_in = zeros(ComplexF64, 4, 3)
    k_in[:, 1] .= ξ

    for (i, q_i) in enumerate(q_in)
        k_in[i, 3] = q_i
    end
    
    k_in ./= c_0
    
    E_forward_in_p =  γ_in[1, :]  # p-polarized incident electric field
    E_forward_in_s =  γ_in[2, :]  # s-polarized incident electric field
    # E_backward_in_p = γ_in[3, :]
    # E_backward_in_s = γ_in[4, :]

    E_out_p = t_coefs[1] * γ_out[1, :] + t_coefs[2] * γ_out[2, :]
    E_out_s = t_coefs[3] * γ_out[1, :] + t_coefs[4] * γ_out[2, :]

    E_ref_p = r_coefs[1] * γ_in[3, :] + r_coefs[2] * γ_in[4, :]
    E_ref_s = r_coefs[3] * γ_in[3, :] + r_coefs[4] * γ_in[4, :]

    S_in_p = real( 0.5 * E_forward_in_p × conj(k_in[1, :] × E_forward_in_p) )
    S_in_s = real( 0.5 * E_forward_in_s × conj(k_in[2, :] × E_forward_in_s) )

    k_out = zeros(ComplexF64, 4, 3)
    k_out[:, 1] .= ξ

    for (i, q_i) in enumerate(q_out)
        k_out[i, 3] = q_i
    end

    k_out ./= c_0

    S_out_p  = real( 0.5 * E_out_p × conj(k_out[1, :] × E_out_p) )
    S_out_s  = real( 0.5 * E_out_s × conj(k_out[2, :] × E_out_s) )
    S_refl_p = real( 0.5 * E_ref_p × conj(k_out[3, :] × E_ref_p) )
    S_refl_s = real( 0.5 * E_ref_s × conj(k_out[4, :] × E_ref_s) )


    return Poynting(S_out_p, S_in_p, S_out_s, S_in_s, S_refl_p, S_refl_s)
end

"""
NOT YET IMPLEMENTED
"""
function electric_field(γ, T, P)
    E_t_po = γ[1, :]
    E_t_se = γ[2, :]
    E_r_po = γ[3, :]

end

"""
Ratio of the absolution squares of two components
used to evaluate if a material is birefringent.
"""
function abs_ratio(a, b)
    return abs2(a) / (abs2(a) + abs2(b))
end

"""
For the four modes (two transmitting and two reflecting), the ratio

        C = |Ex|^2 / (|Ex|^2 + |Ey|^2)

          = |Ψ_1|^2 / (|Ψ_1|^2 + |Ψ_3|^2)

is evaluated. Recall that the values for the electric field are contained
in the eigenvector matrix, Ψ.

If the layer material is birefringent, there will be anisotropy in the
dielectric tensor. If this is the case, the x and y components of the Poynting vector needs to 
be analyzed (eqn 15 in Passler et al., 2017):

        C = |Sx|^2 / (|Sx|^2 + |Sy|^2).

If there is no birefringence, then the electric field is analyzed.
This analysis follows Li et al., 1988.
DOI: 10.1364/AO.27.001334

and the use of the Poynting vector is from Passler et al., 2017, 2019
DOI: 10.1364/JOSAB.34.002128
DOI: 10.1364/JOSAB.36.003246
"""
function evaluate_birefringence(Ψ, S, t_modes, r_modes)

    C_q1 = abs_ratio(S[1, t_modes[1]], S[2, t_modes[1]])
    C_q2 = abs_ratio(S[1, t_modes[2]], S[2, t_modes[2]])
    
    # Note: isapprox(NaN, NaN) = false, which is important
    # in the case that both Sx and Sy are zero.
    if isapprox(C_q1, C_q2)

        if C_q2 > C_q1
            reverse!(t_modes)
        end

        C_q3 = abs_ratio(S[1, r_modes[1]], S[2, r_modes[1]])
        C_q4 = abs_ratio(S[1, r_modes[2]], S[2, r_modes[2]])

        if C_q4 > C_q3
            reverse!(r_modes)
        end
    else
        C_q1 = abs_ratio(Ψ[1, t_modes[1]], Ψ[3, t_modes[1]])
        C_q2 = abs_ratio(Ψ[1, t_modes[2]], Ψ[3, t_modes[2]])

        if C_q2 > C_q1
            reverse!(t_modes)
        end
        C_q3 = abs_ratio(Ψ[1, r_modes[1]], Ψ[3, r_modes[1]])
        C_q4 = abs_ratio(Ψ[1, r_modes[2]], Ψ[3, r_modes[2]])

        if C_q4 > C_q3
            reverse!(r_modes)
        end

    end
    return t_modes, r_modes
end


function calculate_q(Δ, a)

    q_unsorted, Ψ_unsorted = eigen(Δ)

    transmitted_mode = Vector{Integer}(undef, 2)
    reflected_mode = Vector{Integer}(undef, 2)

    kt = 1
    kr = 1

    if isreal(q_unsorted)
        for m in 1:4
            if real(q_unsorted[m]) >= 0.0
                transmitted_mode[kt] = m
                kt += 1
            else
                reflected_mode[kr] = m
                kr += 1
            end
        end
    else
        for m in 1:4
            if imag(q_unsorted[m]) >= 0.0
                transmitted_mode[kt] = m
                kt += 1
            else
                reflected_mode[kr] = m
                kr += 1
            end
        end
    end

    S_unsorted = poynting(Ψ_unsorted, a)
    t, r = evaluate_birefringence(Ψ_unsorted, S_unsorted, transmitted_mode, reflected_mode)

    q = [q_unsorted[t[1]], q_unsorted[t[2]], q_unsorted[r[1]], q_unsorted[r[2]]]
    S = [S_unsorted[t[1]], S_unsorted[t[2]], S_unsorted[r[1]], S_unsorted[r[2]]]

    return q, S
end

"""
These are factors that belong to the electric field calculated
such that singularities are identified and removed.
q[1] and q[2] are forward-traveling modes and
q[3] and q[4] are backward-traveling modes.

This is based on the work in:
Xu et al. Optical degeneracies in anisotropic layered media:
Treatment of singularities in a 4x4 matrix formalism, 2000.
DOI: 10.1103/PhysRevB.61.1740
"""
function calculate_γ(ξ, q, ε, μ)

    γ = Array{ComplexF64}(undef, 4, 3)

    γ[1,1] =  1
    γ[2,2] =  1
    γ[4,2] =  1
    γ[3,1] = -1

    if isapprox(q[1], q[2])

        γ[1,2] = 0
        γ[2,1] = 0
    else
        γ[1,2] = ( μ * ε[2,3] * (μ * ε[3,1] + ξ * q[1])  -  μ * ε[2,1] * (μ * ε[3,3] - ξ^2)) / ( (μ * ε[3,3] - ξ^2) * (μ * ε[2,2] - ξ^2 - q[1]^2) - μ^2 * ε[2,3] * ε[3,2] )
        γ[2,1] = ( μ * ε[3,2] * (μ * ε[1,3] + ξ * q[2])  -  μ * ε[1,2] * (μ * ε[3,3] - ξ^2)) / ( (μ * ε[3,3] - ξ^2) * (μ * ε[1,1] - q[2]^2) - (μ * ε[1,3] + ξ * q[2]) * (μ * ε[3,1] + ξ * q[2]) )
    end

    γ[1,3] = ( -μ * ε[3,1] - ξ * q[1] - μ * ε[3,2] * γ[1,2] ) / (μ * ε[3,3] - ξ^2)
    γ[2,3] = (-(μ * ε[3,1] + ξ * q[2]) * γ[2,1] - μ * ε[3,2]) / (μ * ε[3,3] - ξ^2)

    if isapprox(q[3], q[4])

        γ[3,2] = 0.0
        γ[4,1] = 0.0
    else
        γ[3,2] = ( μ * ε[2,1] * (μ * ε[3,3] - ξ^2) - μ * ε[2,3] * (μ * ε[3,1] + ξ * q[3]) ) / ( (μ * ε[3,3] - ξ^2) * (μ * ε[2,2] - ξ^2 - q[3]^2) - μ^2 * ε[2,3] * ε[3,2] )
        γ[4,1] = ( μ * ε[3,2] * (μ * ε[1,3] + ξ * q[4]) - μ * ε[1,2] * (μ * ε[3,3] - ξ^2) ) / ( (μ * ε[3,3] - ξ^2) * (μ * ε[1,1] - q[4]^2) - (μ * ε[1,3] + ξ * q[4]) * (μ * ε[3,1] + ξ * q[4]) )
    end

    γ[3,3] = ( μ * ε[3,1] + ξ * q[3] + μ * ε[3,2] * γ[3,2] ) / (μ * ε[3,3] - ξ^2)
    γ[4,3] = ( -(μ * ε[3,1] + ξ * q[4]) * γ[4,1] - μ * ε[3,2] ) / (μ * ε[3,3] - ξ^2)

    # normalize γ
    for i in 1:4
        Z = √( γ[i, :] ⋅ γ[i, :]' )
        for j in 1:3
            γ[i, j] /= Z
        end
    end

    return γ
end

function dynamical_matrix(ξ, q, γ, μ)

    A = Matrix{ComplexF64}(undef, 4, 4)

    A[1, :] =  γ[:, 1]
    A[2, :] =  γ[:, 2]
    A[3, :] = (γ[:, 1] .* q .- ξ * γ[:, 3]) ./ μ
    A[4, :] =  γ[:, 2] .* q ./ μ

    return A
end

function propagation_matrix(ω, q, z)
    return Diagonal(exp.(-im * ω * q * z / c_0))
end

function transfermatrix(ω, ξ, q, γ, μ, d)

    A_i = dynamical_matrix(ξ, q, γ, μ)
    P_i = propagation_matrix(ω, q, d)
    T_i = A_i * P_i * inv(A_i)

    return A_i, P_i, T_i
end

function layer_params(ω, ξ, n, μ, d)

    ε_i = dielectric_constant(n)
    ε = dielectric_tensor(ε_i, ε_i, ε_i)

    M = zeros(ComplexF64, 6, 6)
    M[1:3, 1:3] = ε
    M[4:6, 4:6] = Diagonal([μ, μ, μ])
    
    a = construct_a(ξ, M)
    Δ = constructΔ(ξ, M, a)
    q, S = calculate_q(Δ, a)
    γ = calculate_γ(ξ, q, ε, μ)
    A, P, T = transfermatrix(ω, ξ, q, γ, μ, d)

    return A, P, T, γ, q
end

"""
Given a new set of wavelengths, interpolate the 
complex refractive index values for the input Layer
and return a new Layer with the new λ, n, and κ.
The new wavelengths must not extend beyond the 
domain of the existing wavelengths in the Layer (i.e. no extrapolation).
"""
function interp_data(layer::Layer, λs)

    if length(layer.λ) == 1
        n = fill(layer.n[1], length(λs))
        κ = fill(layer.κ[1], length(λs))

        return Layer(layer.material, layer.thickness, λs, n, κ)
    else
        interp_n = LinearInterpolation(layer.n, layer.λ)
        interp_κ = LinearInterpolation(layer.κ, layer.λ)

        n = interp_n.(λs)
        κ = interp_κ.(λs)

        return Layer(layer.material, layer.thickness, λs, n, κ)
    end
end

function calculate_Γ_S(s::Structure, θ::Float64)

    superstrate = s.layers[1]
    substrate = s.layers[end]

    μ = 1.0 + 0.0im

    ε_0in = dielectric_constant(superstrate)
    ξs = @. √(ε_0in) * sin(θ)

    Λ_1324 = [1 0 0 0;
              0 0 1 0;
              0 1 0 0;
              0 0 0 1]

    Γs = Vector{Matrix{ComplexF64}}([])
    Ss = Vector{Poynting}([])
    γ_0s = []
    γ_fs = []
    q_0s = []
    q_fs = []

    ωs = 2π * c_0 ./ s.λ

    for (i, ω) in enumerate(ωs)

        Γ = I

        A_0, P_0, T_0, γ_0, q_0 = layer_params(ω, ξs[i], superstrate.n[i] + superstrate.κ[i] * im, μ, superstrate.thickness)
        A_f, P_f, T_f, γ_f, q_f = layer_params(ω, ξs[i], substrate.n[i] + substrate.κ[i] * im, μ, substrate.thickness)

        for layer in s.layers[2:end - 1]

            n = layer.n[i] + layer.κ[i] * im
            A_i, P_i, T_i, γ_i, q_i = layer_params(ω, ξs[i], n, μ, layer.thickness)
            Γ *= T_i
        end

        Γ = inv(Λ_1324) * inv(A_0) * Γ * A_f * Λ_1324
        r, R, t, T = tr_from_Γ(Γ)
        S = poynting(ξs[i], q_0, q_f, γ_0, γ_f, t, r)

        push!(Γs, Γ)
        push!(Ss, S)
    end

    E = electric_field(γ, T, )

    return Γs, Ss
end

"""
Calculate reflectance and transmittance for the total structure.
This takes the matrix Γ*, but for brevity we call it Γ in this function.

This follows the formalism in:
Yeh, Electromagnetic propagation in birefringent layered media, 1979.
DOI: 10.1364/JOSA.69.000742
"""
function tr_from_Γ(Γ::Matrix{ComplexF64})

    d = Γ[1,1] * Γ[3,3] - Γ[1,3] * Γ[3,1]

    rpp = (Γ[2,1] * Γ[3,3] - Γ[2,3] * Γ[3,1]) / d
    rss = (Γ[1,1] * Γ[4,3] - Γ[4,1] * Γ[1,3]) / d
    rps = (Γ[4,1] * Γ[3,3] - Γ[4,3] * Γ[3,1]) / d
    rsp = (Γ[1,1] * Γ[2,3] - Γ[2,1] * Γ[1,3]) / d

    tpp =  Γ[3,3] / d
    tss =  Γ[1,1] / d
    tps = -Γ[3,1] / d
    tsp = -Γ[1,3] / d

    Rpp = abs2(rpp)
    Rss = abs2(rss)
    Rps = abs2(rps)
    Rsp = abs2(rsp)

    Tpp = abs2(tpp)
    Tss = abs2(tss)
    Tps = abs2(tps)
    Tsp = abs2(tsp)

    r = [rpp, rps, rss, rsp]
    R = [Rpp, Rss, Rsp, Rps]

    t = [tpp, tps, tsp, tss]
    T = [Tpp, Tss, Tsp, Tps]

    return r, R, t, T
end

function tr_from_Γ(Γs::Vector{Matrix{ComplexF64}})

    rs = []
    Rs = []
    ts = []
    Ts = []
    for Γ in Γs
        r, R, t, T = tr_from_Γ(Γ)
        push!(rs, r)
        push!(Rs, R)
        push!(ts, t)
        push!(Ts, T)
    end

    return rs, Rs, ts, Ts
end

"""
Calculate transmittance from the Poynting vector struct,
which contains incident and transmitted energy for both
p-polarized and s-polarized waves.
"""
function tr_from_poynting(S::Poynting)

    Tpp = S.out_p[3] / S.in_p[3]
    Tss = S.out_s[3] / S.in_s[3]

    Rpp = - S.refl_p[3] / S.in_p[3]
    Rss = - S.refl_s[3] / S.in_s[3]

    return Tpp, Tss, Rpp, Rss
end

function tr_from_poynting(Ss::Vector{Poynting})
    Tpps = []
    Tsss = []
    Rpps = []
    Rsss = []
    for S in Ss
        Tpp, Tss, Rpp, Rss = tr_from_poynting(S)
        push!(Tpps, Tpp)
        push!(Tsss, Tss)
        push!(Rpps, Rpp)
        push!(Rsss, Rss)
    end
    return Tpps, Tsss, Rpps, Rsss
end





"""
Iterate through each angle provided in the structure
to find the reflectance and transmittance.
"""
function angle_resolved(s::Structure)

    Rpp_spectrum = Matrix{Float64}(undef, length(s.θ), length(s.λ))
    Rss_spectrum = Matrix{Float64}(undef, length(s.θ), length(s.λ))
    Tpp_spectrum = Matrix{Float64}(undef, length(s.θ), length(s.λ))
    Tss_spectrum = Matrix{Float64}(undef, length(s.θ), length(s.λ))
    Γs = []

    for (i, θ) in enumerate(s.θ)

        Γs, Ss = calculate_Γ_S(s, θ)
        rs, Rs, ts, Ts = tr_from_Γ(Γs)

        #TODO: Reflectivity from the Poynting vector has a bug in it.
        Tpp, Tss, Rpp, Rss = tr_from_poynting(Ss)
        Rpp_spectrum[i, :] = [R[1] for R in Rs]
        Rss_spectrum[i, :] = [R[2] for R in Rs]
        # Rpp_spectrum[i, :] = Rpp
        # Rss_spectrum[i, :] = Rss
        Tpp_spectrum[i, :] = Tpp
        Tss_spectrum[i, :] = Tss
    end

    return Rpp_spectrum, Rss_spectrum, Tpp_spectrum, Tss_spectrum, Γs
end




function transmitted_power()

end

"""
Initialize the structure after all layers have been added.
"""
function initialize!(λs, structure::Structure)
    layers = Layer[]
    for layer in structure.layers
        new_layer = interp_data(layer, λs)
        push!(layers, new_layer)
    end
    return Structure(layers, λs, structure.θ)
end

"""
Print each layer and it's thickness in a somewhat 
visually useful way.
"""
function printstruct(s::Structure)
    print("\n")
    for layer in s.layers

        print("-"^30, "\n")
        print("    $(layer.material), d = $(layer.thickness * 10^6) μm\n")

    end
    print("-"^30, "\n")
end