using Peaks
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

function find_resonance(spectra, atol=1e-3)
    for i in eachindex(spectra)
        pks, vals = findmaxima(spectra[i, :])
        if length(pks) == 2 && isapprox(vals[1], vals[2], atol=atol)
            return i, pks
        end
    end
end

function draw_index_profile(ax, indices, thicknesses)
    prev_x = 0
    prev_n = indices[1]
    for (i, n) in enumerate(indices)
        current_x = sum(thicknesses[1:i])
        lines!(ax, [prev_x, current_x], [n, n], color = :black, linewidth = 0.5)  # Plot the horizontal line
        if i > 1
            lines!(ax, [prev_x, prev_x], [prev_n, n], color = :black, linewidth = 0.5)  # Plot the vertical line
        end
        prev_x = current_x
        prev_n = n
    end
end

n_air = RefractiveMaterial("other", "air", "Ciddor")
n_tio2 = RefractiveMaterial("main", "TiO2", "Kischkat")
n_sio2 = RefractiveMaterial("main", "SiO2", "Kischkat")

λ_0 = 5.0
λs = range(4.8, 5.2, length = 200)
frequencies = 10^4 ./ λs
θs = range(0, 30, length = 50)

# absorbing material
n_bg = 1.4
A_0 = 3000.0
ω_0 = 10^4 / λ_0  # cm^-1
Γ_0 = 5
p0 = [A_0, ω_0, Γ_0]
ε1 = dielectric_real(frequencies, p0) .+ n_bg^2
ε2 = dielectric_imag(frequencies, p0)
n_medium = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) + ε1) / 2)
k_medium = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) - ε1) / 2)

# dbr layers
t_tio2 = λ_0 / (4 * n_tio2(λ_0))
t_sio2 = λ_0 / (4 * n_sio2(λ_0))
t_cav = 1 * λ_0  / n_bg + 0.1  # Slightly offset the cavity length to get negative detuning


air = Layer(n_air, 2.0);
tio2 = Layer(n_tio2, t_tio2);
sio2 = Layer(n_sio2, t_sio2);
absorber = Layer(λs, n_medium, k_medium, t_cav)
absorber.dispersion(4.9)

nperiods = 6
unit = [tio2, sio2]
layers = [air, repeat(unit, nperiods)..., absorber, repeat(reverse(unit), nperiods)..., air];

res = angle_resolved(λs, deg2rad.(θs), layers)
θ_idx, peaks = find_resonance(res.Tpp)
T_plot = res.Tpp[θ_idx, :]
θ_plot = θs[θ_idx]

field1 = electric_field(λs[peaks[1]], layers)
field2 = electric_field(λs[peaks[2]], layers)

##

fig = Figure(size = (750, 600))
ax1 = Axis(fig[1, 1], title = "Polariton dispersion",
    xlabel = "Incidence angle (°)",
    ylabel = "Frequency (cm⁻¹)")
heatmap!(θs, νs, res.Tpp, colormap = :deep)


ax2 = Axis(fig[2, 1], title = "Normal mode splitting at (θ ≈ $(round(Int, θ_plot)))°",
    xlabel = "Frequency (cm⁻¹)",
    ylabel = "Transmittance")
lines!(frequencies, T_plot)
scatter!(frequencies[peaks], T_plot[peaks], color = :red, marker = 'x', markersize = 15)

ax3 = Axis(fig[1, 2], ylabel = "Photon field")
lines!(field1.z .* 1e3, real(field1.p[1, :]), label = "LP")
lines!(field2.z .* 1e3, real(field2.p[1, :]), label = "UP")
hlines!(0, color = :black, linestyle = :dash, linewidth = 0.5)

n1 = round(real(sio2.dispersion(λ_0)), digits=2)
n2 = round(real(tio2.dispersion(λ_0)), digits=2)

ax4 = Axis(fig[2, 2], xlabel = "Distance (μm)", ylabel = "Refractive index", yticks = [n1, n2])

# Draw DBR cavity structure
refractive_indices = [real(layer.dispersion(λ_0)) for layer in layers[2:end-1]]
thicknesses = [layer.thickness for layer in layers[2:end-1]]

draw_index_profile(ax4, refractive_indices, thicknesses)

axislegend(ax3)
hidexdecorations!(ax3)
hideydecorations!(ax3, label = false, ticks = false, ticklabels = false)
hidexdecorations!(ax4, label = false, ticks = false, ticklabels = false)
hideydecorations!(ax4, label = false, ticks = false, ticklabels = false)

fig