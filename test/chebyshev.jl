using Test
import StarAlgebras as SA
using MultivariateBases
import MultivariateBases as MB
using DynamicPolynomials

@testset "StarAlgebras" begin
    @polyvar x
    full = MB.FullBasis{MB.Chebyshev}(x)
    sub = SA.SubBasis(full, [[1]])
    a = MB.algebra_element(MB.sparse_coefficients(1 // 1 * x), full)
    @test a == MB.convert_basis(full, MB.algebra_element(x))
    b = a * a
    expected =
        MB.algebra_element(MB.sparse_coefficients((1 // 2) * (1 + x^2)), full)
    @test b == expected
    c = b * b
    @test SA.coeffs(c) ==
          MB.sparse_coefficients(3 // 8 + 1 // 2 * x^2 + 1 // 8 * x^4)
    c1 = MB.convert_basis(full, MB.algebra_element(constant_monomial(x)))
    c2 = MB.convert_basis(full, MB.algebra_element(x^0))
    @test a * c1 == a * c2
end

@testset "Orthogonal" begin
    orthogonal_test(
        MB.Chebyshev,
        x -> [1, x, 2x^2 - 1, 4x^3 - 3x, 8x^4 - 8x^2 + 1],
        true,
    )
    univ_orthogonal_test(
        MB.Chebyshev,
        i -> π * (1 / 2 + (0^i) / 2);
        atol = 1e-12,
    )
    orthogonal_test(
        MB.ChebyshevSecondKind,
        x -> [1, 2x, 4x^2 - 1, 8x^3 - 4x, 16x^4 - 12x^2 + 1],
        true,
    )

    univ_orthogonal_test(MB.ChebyshevSecondKind, i -> π / 2; atol = 1e-12)
end

@testset "API degree = $degree" for degree in 0:3
    api_test(MB.Chebyshev, degree)
    api_test(MB.ChebyshevSecondKind, degree)
end

@testset "Coefficients" begin
    coefficient_test(
        MB.Chebyshev,
        reverse([
            0.0625,
            0.0625,
            0.0625,
            -0.25,
            0.0625,
            -0.3125,
            -0.3125,
            0.625,
        ]),
    )
    coefficient_test(
        MB.ChebyshevSecondKind,
        reverse([
            0.015625,
            0.015625,
            0.015625,
            -0.09375,
            0.015625,
            -0.109375,
            -0.109375,
            0.875,
        ]),
    )
end
