using Test
using MultivariateBases
using DynamicPolynomials

@testset "StarAlgebras" begin
    @polyvar x
    full = MB.FullBasis{MB.Trigonometric}(x)
    a = MB.algebra_element(MB.sparse_coefficients(1//1 * x), full)
    b = a * a
    @test b.coeffs == MB.sparse_coefficients(1 // 2 + 1 // 2 * x^3)
    c = b * b
    @test c.coeffs ==
          MB.sparse_coefficients(3 // 8 + 1 // 2 * x^3 + 1 // 8 * x^7)
    d = MB.algebra_element(
        MB.sparse_coefficients(1//1 * constant_monomial(typeof(x))),
        full,
    )
    e = MB.algebra_element(MB.sparse_coefficients(1//1 * x^0), full)
    @test a * d == a * e
end
