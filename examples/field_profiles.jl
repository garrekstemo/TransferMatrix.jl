# Field profiles (E and H): three things the magnetic field shows that |E| cannot.
#
# 1. TM duality: for p-polarized incidence the magnetic field is a single transverse
#    component (Hy) while the electric field genuinely splits into Ex and Ez.
# 2. E/H complementarity: in a resonant cavity the |E| antinodes coincide with |H| nodes;
#    the H antinodes sit at the mirrors where E vanishes (relevant for magnetic-dipole
#    coupling — a molecule at the E maximum sits in a magnetic null).
# 3. Sheet H-jump: across a 2D conductive Sheet the in-plane H jumps by the surface
#    current σ̃·E while tangential E stays continuous — an H-only observable.
#
# hfield returns the impedance-normalized magnetic field H̃ = Z₀·H_SI, so |E| and |H̃|
# overlay on one scale. Run from the examples environment:
#   julia --project=examples examples/field_profiles.jl

using TransferMatrix
using RefractiveIndex
using CairoMakie

# Panel A: TM duality (oblique p-incidence)
airA  = Layer(λ -> complex(1.0), 0.3)
filmA = Layer(λ -> complex(2.0), 0.3)
subA  = Layer(λ -> complex(1.5), 0.0)
layersA = [airA, filmA, subA]
θ = deg2rad(45); λA = 0.6; dz = 0.002
efA = efield(λA, layersA; θ=θ, dz=dz)
hfA = hfield(λA, layersA; θ=θ, dz=dz)
ExA = abs.(efA.p[1, :]); EzA = abs.(efA.p[3, :])
HxA = abs.(hfA.p[1, :]); HyA = abs.(hfA.p[2, :]); HzA = abs.(hfA.p[3, :])
mA = maximum(vcat(ExA, EzA, HyA))
println("Panel A  |Ex|=", round(maximum(ExA),digits=3), " |Ez|=", round(maximum(EzA),digits=3),
        " |H̃y|=", round(maximum(HyA),digits=3), " | |H̃x|,|H̃z|=", round(maximum(HxA),digits=3), ",", round(maximum(HzA),digits=3))

# Panel B: E/H complementarity in an IR Au/air/Au cavity at resonance
au    = RefractiveMaterial("main", "Au", "Rakic-LD")
air_m = RefractiveMaterial("other", "air", "Ciddor")
λ0 = 5.0; t_middle = λ0 / 2
air = Layer(air_m, t_middle); auL = Layer(au, 0.01)
layers = [air, auL, air, auL, air]
λs = range(2.0, 6.0, length = 500)
Tpp = [transfer(λ, layers).Tpp for λ in λs]
mask = (λs .>= 4.9) .& (λs .<= 5.5)
λres = collect(λs)[mask][argmax(Tpp[mask])]
ef = efield(λres, layers; dz=0.002)
hf = hfield(λres, layers; dz=0.002)
Emag = vec(sqrt.(sum(abs2, ef.p, dims=1))); Emag ./= maximum(Emag)
Hmag = vec(sqrt.(sum(abs2, hf.p, dims=1))); Hmag ./= maximum(Hmag)
iE = argmax(Emag)
println("Panel B  resonance λ=", round(λres,digits=4), " μm; E antinode at z=", round(ef.z[iE],digits=3),
        " μm has |H|=", round(Hmag[iE],digits=3))

# Panel C: 2D Sheet H-jump (s-incidence)
n0 = 1.0
airC = Layer(λ -> complex(n0), 0.0)
spacerL = Layer(λ -> complex(1.0), 0.5)
spacerR = Layer(λ -> complex(1.0), 0.5)
subC = Layer(λ -> complex(n0), 0.0)
layersC = [airC, spacerL, spacerR, subC]
σ_s = 1.0e-3 + 0.0im
sheetsC = Dict(2 => Sheet(σ_s))
efC = efield(0.6, layersC; dz=5e-4, sheets=sheetsC)
hfC = hfield(0.6, layersC; dz=5e-4, sheets=sheetsC)
EyC = abs.(efC.s[2, :]); HxC = abs.(hfC.s[1, :])
zint = efC.boundaries[2]
jb = findlast(z -> z ≤ zint, efC.z); ja = jb + 1
println("Panel C  ΔH̃x=", round(abs(hfC.s[1,ja]-hfC.s[1,jb]),digits=4),
        " (Z₀σ̃·Ey=", round(abs(376.730313668*σ_s*efC.s[2,jb]),digits=4),
        "); |Ey| jump=", round(abs(EyC[ja]-EyC[jb]),digits=5))

fig = Figure(size = (820, 940))

axA = Axis(fig[1, 1]; xlabel = "z (μm)", ylabel = "|field| (norm.)",
    title = "TM (p) incidence: H is single-component (H̃y); E is two-component (Ex, Ez)")
vlines!(axA, efA.boundaries; color = :gray, linestyle = :dot, linewidth = 1)
lines!(axA, efA.z, ExA ./ mA; color = :dodgerblue3, linewidth = 2, label = "|Ex|")
lines!(axA, efA.z, EzA ./ mA; color = :dodgerblue1, linewidth = 2, linestyle = :dash, label = "|Ez|")
lines!(axA, hfA.z, HyA ./ mA; color = :firebrick3, linewidth = 2.5, label = "|H̃y|")
lines!(axA, hfA.z, HxA ./ mA; color = (:gray, 0.7), linewidth = 1, label = "|H̃x|,|H̃z| ≈ 0")
lines!(axA, hfA.z, HzA ./ mA; color = (:gray, 0.7), linewidth = 1)
axislegend(axA; position = :lt, framevisible = true)

axB = Axis(fig[2, 1]; xlabel = "z (μm)", ylabel = "|field| (each norm.)",
    title = "Au/air/Au IR cavity at resonance (λ=$(round(λres,digits=2)) μm): |E| antinode = |H| node")
vlines!(axB, ef.boundaries; color = :gray, linestyle = :dot, linewidth = 1)
lines!(axB, ef.z, Emag; color = :dodgerblue3, linewidth = 2.5, label = "|E|")
lines!(axB, hf.z, Hmag; color = :firebrick3, linewidth = 2.5, label = "|H̃| (= Z₀|H_SI|)")
text!(axB, ef.z[iE], 1.0; text = "E antinode\n(H null here)", align = (:center, :top), fontsize = 10, color = :dodgerblue4)
axislegend(axB; position = :lt, framevisible = true)

axC = Axis(fig[3, 1]; xlabel = "z (μm)", ylabel = "|field| (norm.)",
    title = "2D conductive Sheet (s-pol): tangential E continuous, H̃ jumps by σ̃·E")
mC = maximum(vcat(EyC, HxC))
vlines!(axC, [zint]; color = :seagreen, linestyle = :dash, linewidth = 1.5)
lines!(axC, efC.z, EyC ./ mC; color = :dodgerblue3, linewidth = 2.5, label = "|Ey| (continuous)")
lines!(axC, hfC.z, HxC ./ mC; color = :firebrick3, linewidth = 2.5, label = "|H̃x| (jumps)")
text!(axC, zint, 0.05; text = "sheet", align = (:center, :bottom), fontsize = 10, color = :seagreen)
axislegend(axC; position = :lt, framevisible = true)

outpath = joinpath(@__DIR__, "..", "docs", "src", "assets", "examples", "field_profiles.png")
mkpath(dirname(outpath))
save(outpath, fig)
println("saved ", outpath)
