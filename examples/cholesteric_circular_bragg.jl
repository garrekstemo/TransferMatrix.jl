# Cholesteric selective reflection — the textbook showcase for the circular basis.
#
# A cholesteric (chiral nematic) liquid crystal is a stack of birefringent layers whose
# optic axis rotates helically about the propagation (z) axis with a fixed pitch P (the
# distance for a full 360° turn). It is NOT intrinsically optically active — each slice is
# an ordinary uniaxial crystal — yet the helical STRUCTURE makes it reflect circularly
# polarized light of one handedness almost totally inside a Bragg band, while passing the
# opposite handedness. This is "circular Bragg" / selective reflection, the effect behind
# reflective cholesteric-LC polarizers and color filters.
#
#   de Gennes & Prost, The Physics of Liquid Crystals (1993), §6.3;
#   D. W. Berreman & T. J. Scheffer, Phys. Rev. Lett. 25, 577 (1970).
#
# Standard reference you can check by hand (normal incidence):
#   band center     λ0 = n̄ P            with n̄ = (n_o + n_e)/2
#   band edges      n_o P  <  λ  <  n_e P
#   band width      Δλ = Δn P           with Δn = n_e − n_o
# Inside the band the co-handed reflectance → 1 for a thick film; the cross-handed
# reflectance ≈ 0. The LINEAR basis hides this: Rpp = Rss ≈ 0.5 in the band (half of any
# linear/unpolarized beam is reflected), revealing nothing about handedness. The circular
# basis splits it cleanly — that is why `basis=:circular` exists.
#
# We have no cholesteric primitive yet (issue #87/#94), so we build the helix by hand as a
# fine stack of rotated uniaxial `Layer`s — the same piecewise-uniform approximation those
# issues will wrap. The slices are wavelength-independent, so we build the stack ONCE and
# reuse it across the spectrum.
#
# Run from the examples environment (already dev-pointed at this checkout):
#   julia --project=examples examples/cholesteric_circular_bragg.jl

using TransferMatrix
using CairoMakie

# --- cholesteric parameters (a typical visible-band chiral LC) ---
n_o = 1.5                       # ordinary index (⊥ optic axis)
n_e = 1.7                       # extraordinary index (∥ optic axis)
P = 0.40                        # helical pitch (μm), full 360° turn
n_bar = (n_o + n_e) / 2         # mean index
Δn = n_e - n_o

n_pitch = 30                    # number of full helical turns (film ≈ 12 μm thick)
slices_per_pitch = 24           # angular resolution of the helix (thin slices)

# Build the helix: a uniaxial slab (optic axis along z in its own frame) laid into the
# xy-plane with euler=(φ, π/2, 0) — the ZYZ angle θ=π/2 tilts the optic axis into the
# plane and φ sets its azimuth. Stepping φ = q·z (q = 2π/P) winds the director through z.
# Sign of q = helix handedness.
function cholesteric(q)
    d_slice = P / slices_per_pitch
    n_slices = n_pitch * slices_per_pitch
    slices = [Layer(λ -> n_o + 0im, λ -> n_o + 0im, λ -> n_e + 0im, d_slice;
                    euler = (q * (i - 0.5) * d_slice, π/2, 0.0)) for i in 1:n_slices]
    # Ambient index-matched to n̄ on both sides: isolates the chiral Bragg effect from
    # ordinary Fresnel boundary reflections, so the spectrum is a clean comparison to the
    # analytic band. (Isotropic ambient also keeps us clear of the anisotropic-ambient
    # limitation, issue #71.)
    ambient = Layer(λ -> n_bar + 0im, 0.0)
    return [ambient; slices; ambient]
end

q = 2π / P
helix = cholesteric(+q)         # one handedness

# --- spectra: circular co/cross reflectance, plus linear for contrast ---
λs = range(0.50, 0.80, length = 300)
R_LL = similar(collect(λs))     # left-circular reflected from left-circular incident
R_RR = similar(collect(λs))     # right-circular reflected from right-circular incident
R_lin = similar(collect(λs))    # TOTAL linear reflectance for p-input: R_pp + R_sp
for (i, λ) in enumerate(λs)
    rc = transfer(λ, helix; basis = :circular)
    R_LL[i] = rc.Rll
    R_RR[i] = rc.Rrr
    rl = transfer(λ, helix)
    R_lin[i] = rl.Rpp + rl.Rsp   # reflected co-handed light splits across p and s
end

# --- compare to the standard reference ---
λ0 = n_bar * P
edge_lo, edge_hi = n_o * P, n_e * P
inband = findall(>(0.5), R_LL)                 # measured band: where co-handed R > 0.5
meas_lo, meas_hi = λs[first(inband)], λs[last(inband)]
meas_center = (meas_lo + meas_hi) / 2
icenter = argmin(abs.(λs .- λ0))

println("Cholesteric circular Bragg — TMM vs standard reference")
println("  film: $(n_pitch) pitches × $(slices_per_pitch) slices = $(length(helix)-2) layers, ",
        round(n_pitch * P, digits=2), " μm thick")
println("  n_o=$(n_o)  n_e=$(n_e)  P=$(P) μm")
println()
println("  band CENTER   analytic n̄P          = ", round(λ0, digits=4), " μm")
println("                measured (mid R_LL>0.5) = ", round(meas_center, digits=4), " μm")
println("  band EDGES    analytic [n_oP, n_eP]   = [", round(edge_lo,digits=4), ", ", round(edge_hi,digits=4), "] μm")
println("                measured (R_LL>0.5)     = [", round(meas_lo,digits=4), ", ", round(meas_hi,digits=4), "] μm")
println("  band WIDTH    analytic ΔnP            = ", round(Δn*P, digits=4), " μm   measured ≈ ", round(meas_hi-meas_lo, digits=4), " μm")
println("                (measured > ΔnP because finite-thickness edge sidelobes also exceed 0.5)")
println()
println("  at band center λ0=$(round(λ0,digits=3)) μm:")
println("    R_LL (co-handed)        = ", round(R_LL[icenter], digits=4), "   → reflected")
println("    R_RR (cross-handed)     = ", round(R_RR[icenter], digits=4), "   → transmitted")
println("    R_pp+R_sp (linear, p-in)= ", round(R_lin[icenter], digits=4), "   → linear loses ~half, but cannot reveal the handedness")
println()

# --- the selective-reflection RULE: reverse the helix → reverse the reflected handedness
rc_L = transfer(λ0, cholesteric(+q); basis = :circular)
rc_R = transfer(λ0, cholesteric(-q); basis = :circular)
println("  reflected handedness tracks the helix (selective-reflection rule):")
println("    +q helix:  R_LL=", round(rc_L.Rll, digits=3), "  R_RR=", round(rc_L.Rrr, digits=3))
println("    −q helix:  R_LL=", round(rc_R.Rll, digits=3), "  R_RR=", round(rc_R.Rrr, digits=3))

# --- plot the comparison ---
fig = Figure(size = (760, 460))
ax = Axis(fig[1, 1];
    xlabel = "Wavelength (μm)", ylabel = "Reflectance",
    title = "Cholesteric selective reflection (n_o=$(n_o), n_e=$(n_e), P=$(P) μm)")

# analytic Bragg band (the standard reference): shaded span + center line
vspan!(ax, edge_lo, edge_hi; color = (:gray, 0.15))
vlines!(ax, [λ0]; color = :gray, linestyle = :dot, linewidth = 1)

lines!(ax, λs, R_LL; color = :dodgerblue3, linewidth = 2.5, label = "R_LL  (co-handed, circular)")
lines!(ax, λs, R_RR; color = :firebrick3, linewidth = 2.5, label = "R_RR  (cross-handed, circular)")
lines!(ax, λs, R_lin; color = :black, linestyle = :dash, linewidth = 1.5, label = "R_pp+R_sp  (linear total, p-in)")

text!(ax, λ0, 1.04; text = "analytic band\n[n_oP, n_eP]", align = (:center, :bottom),
      fontsize = 11, color = :gray30)
ylims!(ax, -0.03, 1.18)
axislegend(ax; position = :rc, framevisible = true)

outpath = joinpath(@__DIR__, "cholesteric_circular_bragg.png")
save(outpath, fig)
println("\nSaved figure to ", outpath)
