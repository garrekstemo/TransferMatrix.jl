using Test
using TransferMatrix
using RefractiveIndex

@testset "RefractiveIndex extension" begin

    @testset "extension is loaded" begin
        @test Base.get_extension(TransferMatrix, :RefractiveIndexExt) !== nothing
    end

    # Sheet(::RefractiveMaterial, d) and _index_fn(::RefractiveMaterial) are
    # relocated to the extension. The Layer / refractive_index material paths are
    # covered in test/types.jl and test/layer.jl; the sheet paths are exercised here.
    @testset "Sheet from RefractiveMaterial" begin
        sio2 = RefractiveMaterial("main", "SiO2", "Malitson")
        d = 6.5e-4   # μm
        λ = 1.0
        n = refractive_index(sio2)(λ)

        sheet = Sheet(sio2, d)
        σ_expected = -im * (2π * TransferMatrix.c_0 / λ) * TransferMatrix.ε_0 * d * (n^2 - 1)
        @test sheet.conductivity(λ)[1, 1] ≈ σ_expected
        @test sheet.conductivity(λ)[2, 2] ≈ σ_expected
        @test sheet.conductivity(λ)[1, 2] == 0

        # Matches the closure-based constructor for the same index function.
        sheet_fn = Sheet(refractive_index(sio2), d)
        @test sheet.conductivity(λ) ≈ sheet_fn.conductivity(λ)

        # In-plane anisotropic Sheet(nx, ny, d) accepts materials too (via _index_fn).
        au = RefractiveMaterial("main", "Au", "Rakic-LD")
        sheet_aniso = Sheet(sio2, au, d)
        nx = refractive_index(sio2)(λ)
        ny = refractive_index(au)(λ)
        σx = -im * (2π * TransferMatrix.c_0 / λ) * TransferMatrix.ε_0 * d * (nx^2 - 1)
        σy = -im * (2π * TransferMatrix.c_0 / λ) * TransferMatrix.ε_0 * d * (ny^2 - 1)
        @test sheet_aniso.conductivity(λ)[1, 1] ≈ σx
        @test sheet_aniso.conductivity(λ)[2, 2] ≈ σy
    end
end
