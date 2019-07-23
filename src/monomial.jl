"""
    struct MonomialBasis{MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
        monomials::MV
    end

Monomial basis with the monomials of the vector `monomials`.
For instance, `MonomialBasis([1, x, y, x^2, x*y, y^2])` is the monomial basis
for the subspace of quadratic polynomials in the variables `x`, `y`.
"""
struct MonomialBasis{MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
    monomials::MV
end
MonomialBasis(monomials) = MonomialBasis(monovec(monomials))

MultivariatePolynomials.polynomialtype(mb::MonomialBasis{MT}, T::Type) where MT = MultivariatePolynomials.polynomialtype(MT, T)
MultivariatePolynomials.polynomial(f::Function, mb::MonomialBasis) = polynomial(f, mb.monomials)
function MultivariatePolynomials.coefficients(p, ::Type{<:MonomialBasis})
    return coefficients(p)
end
