using RefractiveIndex
using TransferMatrix
using GLMakie

function dielectric_real(ω, p)
    A, ω_0, Γ = p
    return @. A * (ω_0^2 - ω^2) / ((ω^2 - ω_0^2)^2 + (Γ * ω)^2)
end
function dielectric_imag(ω, p)
    A, ω_0, Γ = p
    return @. A * Γ * ω / ((ω^2 - ω_0^2)^2 + (Γ * ω)^2)
end


n_air = RefractiveMaterial("other", "air", "Ciddor")
n_tio2 = RefractiveMaterial("main", "TiO2", "Kischkat")
n_sio2 = RefractiveMaterial("main", "SiO2", "Kischkat")

λ_0 = 5.0
λs = range(4.8, 5.2, length = 300)
νs = 10^4 ./ λs
thicknesses = 0.1:0.01:5

# absorbing material
n_bg = 1.4
A_0 = 2000.0
ω_0 = 10^4 / λ_0
Γ_0 = 5
p0 = [A_0, ω_0, Γ_0]
ε1 = dielectric_real(νs, p0) .+ n_bg^2
ε2 = dielectric_imag(νs, p0)
n_medium = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) + ε1) / 2)
k_medium = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) - ε1) / 2)

# dbr layers
t_tio2 = λ_0 / (4 * n_tio2(λ_0))
t_sio2 = λ_0 / (4 * n_sio2(λ_0))
t_cav = 1 * λ_0  / n_bg + 0.1  # Slightly offset the cavity length to get negative detuning


air = Layer(n_air, 2.0);
tio2 = Layer(n_tio2, t_tio2);
sio2 = Layer(n_sio2, t_sio2);
absorber = Layer(λs, n_medium, k_medium, t_cav);

nperiods = 6
unit = [tio2, sio2]
layers = [air, repeat(unit, nperiods)..., absorber, repeat(reverse(unit), nperiods)..., air];

##
res = tune_thickness(λs, thicknesses, layers, 14)


fig = Figure()
ax = Axis(fig[1, 1], title = "Polariton dispersion",
    xlabel = "Thickness (μm)",
    ylabel = "Frequency (cm⁻¹)")
heatmap!(thicknesses, νs, res.Tpp, colormap = :deep)
fig