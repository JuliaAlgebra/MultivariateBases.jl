struct PlusMinusOne # Reduce modulo `x^2 = 1`
end
function Base.rem(m::AbstractMonomial, ::PlusMinusOne)
    return prod(v^mod(e, 2) for (v, e) in powers(m))
end
function Base.rem(t::AbstractTerm, d::PlusMinusOne)
    return coefficient(t) * rem(monomial(t), d)
end
function Base.rem(p::AbstractPolynomial, d::PlusMinusOne)
    return sum(rem(t, d) for t in terms(p))
end

@testset "PlusMinusOne" begin
    @polyvar x y
    basis =
        MB.QuotientBasis(MB.SubBasis{MB.Monomial}([1, y, x]), PlusMinusOne())
    @test length(basis) == 3
    p = x^3 - 2x^2 * y + 3x^2
    coeffs = [3, -2, 1]
    @test coefficients(p, basis) == coeffs
    a = MB.algebra_element(p)
    @test SA.coeffs(a, basis) == coeffs
    explicit = MB.explicit_basis(a)
    exp = MB.algebra_element(SA.coeffs(a, explicit), explicit)
    @test SA.coeffs(exp, basis) == coeffs
end
