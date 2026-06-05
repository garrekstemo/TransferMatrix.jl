# References

All references used for writing TransferMatrix.jl.

## Transfer matrix method

[1] N. C. Passler, M. Jeannin, and A. Paarmann, Layer-Resolved Absorption of Light in Arbitrarily Anisotropic Heterostructures, Phys. Rev. B 101, 165425 (2020). \
[2] B. Garibello, N. Avilán, J. A. Galvis, and C. A. Herreño-Fierro, On the Singularity of the Yeh 4 × 4 Transfer Matrix Formalism, Journal of Modern Optics 67, 832 (2020). \
[3] N. C. Passler and A. Paarmann, Generalized 4 × 4 Matrix Formalism for Light Propagation in Anisotropic Stratified Media: Study of Surface Phonon Polaritons in Polar Dielectric Heterostructures: Erratum, J. Opt. Soc. Am. B 36, 3246 (2019). \
[4] N. C. Passler and A. Paarmann, Generalized 4 × 4 Matrix Formalism for Light Propagation in Anisotropic Stratified Media: Study of Surface Phonon Polaritons in Polar Dielectric Heterostructures, J. Opt. Soc. Am. B 34, 2128 (2017). \
[5] P. Yeh, Optical Waves in Layered Media (Wiley, 2005). \
[6] W. Xu, L. T. Wood, and T. D. Golding, Optical Degeneracies in Anisotropic Layered Media: Treatment of Singularities in a 4×4 Matrix Formalism, Phys. Rev. B 61, 1740 (2000). \
[7] Z.-M. Li, B. T. Sullivan, and R. R. Parsons, Use of the 4 × 4 Matrix Method in the Optics of Multilayer Magnetooptic Recording Media, Appl. Opt., AO 27, 1334 (1988). \
[8] P. Yeh, Electromagnetic Propagation in Birefringent Layered Media, J. Opt. Soc. Am. 69, 742 (1979). \
[9] D. W. Berreman, Optics in Stratified and Anisotropic Media: 4×4-Matrix Formulation, J. Opt. Soc. Am. 62, 502 (1972).

## 2D conductive sheets (TMDC monolayers)

References supporting the modeling of a 2D material / TMDC monolayer as a surface
optical conductivity sheet (the `Sheet` type, `sheet_matrix`, and the
refractive-index ↔ sheet-conductivity conversion). The oblique-incidence,
anisotropic 4×4 interface matrix used in the code follows the Berreman/Passler
formalism above; the surface-conductivity boundary condition and the
`σ_s = -iωε₀d(n²-1)` conversion follow [10] and [11].

[10] Y. Li and T. F. Heinz, Two-Dimensional Models for the Optical Response of Thin Films, 2D Mater. 5, 025021 (2018). \
[11] B. Majerus, E. Dremetsika, M. Lobet, L. Henrard, and P. Kockaert, Electrodynamics of Two-Dimensional Materials: Role of Anisotropy, Phys. Rev. B 98, 125419 (2018). \
[12] Y. V. Morozov and M. Kuno, Optical Constants and Dynamic Conductivities of Single Layer MoS₂, MoSe₂, and WSe₂, Appl. Phys. Lett. 107, 083103 (2015). \
[13] S. Dufferwiel, S. Schwarz, F. Withers, A. A. P. Trichet, F. Li, M. Sich, O. Del Pozo-Zamudio, C. Clark, A. Nalitov, D. D. Solnyshkov, G. Malpuech, K. S. Novoselov, J. M. Smith, M. S. Skolnick, D. N. Krizhanovskii, and A. I. Tartakovskii, Exciton-Polaritons in van der Waals Heterostructures Embedded in Tunable Microcavities, Nat. Commun. 6, 8579 (2015). \
[14] Y. Li, A. Chernikov, X. Zhang, A. Rigosi, H. M. Hill, A. M. van der Zande, D. A. Chenet, E.-M. Shih, J. Hone, and T. F. Heinz, Measurement of the Optical Dielectric Function of Monolayer Transition-Metal Dichalcogenides: MoS₂, MoSe₂, WS₂, and WSe₂, Phys. Rev. B 90, 205422 (2014). \
[15] R. R. Nair, P. Blake, A. N. Grigorenko, K. S. Novoselov, T. J. Booth, T. Stauber, N. M. R. Peres, and A. K. Geim, Fine Structure Constant Defines Visual Transparency of Graphene, Science 320, 1308 (2008). \
[16] L. A. Falkovsky, Optical Properties of Graphene, J. Phys.: Conf. Ser. 129, 012004 (2008).
