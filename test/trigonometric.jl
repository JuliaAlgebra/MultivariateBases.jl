using Test
using MultivariateBases
import MultivariateBases as MB
import StarAlgebras as SA
using DynamicPolynomials

@testset "StarAlgebras" begin
    @polyvar x
    full = MB.FullBasis{MB.Trigonometric}(x)
    a = MB.algebra_element(SA.coeffs(1//1 * x), full)
    b = a * a
    @test b.coeffs == SA.coeffs(1 // 2 + 1 // 2 * x^3)
    c = b * b
    @test c.coeffs == SA.coeffs(3 // 8 + 1 // 2 * x^3 + 1 // 8 * x^7)
    d = MB.algebra_element(SA.coeffs(1//1 * constant_monomial(x)), full)
    e = MB.algebra_element(SA.coeffs(1//1 * x^0), full)
    @test a * d == a * e
    # sin(ωt) * cos(2ωt) = (sin(3ωt) - sin(ωt)) / 2
    # Exponents: sin(ωt) -> x^2, sin(3ωt) -> x^6, cos(2ωt) -> x^3
    s1 = MB.algebra_element(SA.coeffs(1//1 * x^2), full)
    c2 = MB.algebra_element(SA.coeffs(1//1 * x^3), full)
    @test (s1 * c2).coeffs == SA.coeffs(-1//2 * x^2 + 1//2 * x^6)
end
