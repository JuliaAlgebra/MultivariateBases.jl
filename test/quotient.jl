using Test
using DynamicPolynomials
import MultivariateBases as MB
import StarAlgebras as SA

struct ImRem # Reduce modulo `x^2 = -1`
end
function Base.rem(m::AbstractMonomial, ::ImRem)
    return prod((-1)^rem(e, 2) * v^mod(e, 2) for (v, e) in powers(m))
end
function Base.rem(t::AbstractTerm, d::ImRem)
    return coefficient(t) * rem(monomial(t), d)
end
function Base.rem(p::AbstractPolynomial, d::ImRem)
    return sum(rem(t, d) for t in terms(p))
end

@testset "ImRem" begin
    @polyvar x y
    basis = MB.QuotientBasis(MB.SubBasis{MB.Monomial}([1, y, x]), ImRem())
    @test length(basis) == 3
    p = x^3 - 2x^2 * y + 3x^2
    coeffs = [3, 2, -1]
    @test coefficients(p, basis) == coeffs
    a = MB.algebra_element(p)
    @test SA.coeffs(a, basis) == coeffs
    explicit = MB.explicit_basis(a)
    exp = MB.algebra_element(SA.coeffs(a, explicit), explicit)
    @test SA.coeffs(exp, basis) == coeffs
    q = MB.algebra_element(coeffs, basis)
    @test SA.adjoint_coeffs(q, explicit) == [3, -2, 1]
end

@testset "LinearQuotientBasis" begin
    @polyvar x y
    sub = MB.SubBasis{MB.Monomial}([1, y, x])
    @test_throws DimensionMismatch MB.LinearQuotientBasis(sub, ones(2, 2))
    # Quotient by the span of `x - y`
    basis = MB.LinearQuotientBasis(sub, [1 0; 0 1; 0 1])
    @test length(basis) == 2
    @test sprint(show, basis) == "Quotient of dimension 2 of the span of $sub"
    @test coefficients(2 + 3x + y, basis) == [2, 4]
    @test iszero(coefficients(x - y, basis))
    @test_throws ArgumentError coefficients(x^2, basis)
    a = MB.algebra_element(2 + 3x + y)
    @test SA.coeffs(a, basis) == [2, 4]
    explicit = MB.explicit_basis(a)
    exp = MB.algebra_element(SA.coeffs(a, explicit), explicit)
    @test SA.coeffs(exp, basis) == [2, 4]
    @test SA.adjoint_coeffs(
        [5, 7],
        MB.SubBasis{MB.Monomial}([y, x^2]),
        basis,
    ) == [7, 0]
    q = MB.algebra_element([5, 7], basis)
    @test SA.coeffs(q) == [5, 7]
end
