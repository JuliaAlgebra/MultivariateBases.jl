using Test
using MultivariateBases
using DynamicPolynomials

@testset "StarAlgebras" begin
    @polyvar x
    a = MB.Polynomial{MB.Trigonometric}(x)
    b = a * a
    @test b.coeffs == MB.sparse_coefficients(1 // 2 + 1 // 2 * x^3)
    c = b * b
    @test c.coeffs ==
          MB.sparse_coefficients(3 // 8 + 1 // 2 * x^3 + 1 // 8 * x^7)
    @test a * MB.Polynomial{MB.Trigonometric}(constant_monomial(typeof(x))) ==
          a * MB.Polynomial{MB.Trigonometric}(x^0)
end
