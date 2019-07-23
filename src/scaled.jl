"""
    struct ScaledMonomialBasis{MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
        monomials::MV
    end

*Scaled monomial basis* (see [Section 3.1.5, BPT12]) with the monomials of the vector `monomials`.
Given a monomial ``x^\\alpha = x_1^{\\alpha_1} \\cdots x_n^{\\alpha_n}`` of degree ``d = \\sum_{i=1}^n \\alpha_i``,
the corresponding polynomial of the basis is
```math
{d \\choose \\alpha}^{\\frac{1}{2}} x^{\\alpha} \\quad \\text{ where } \\quad
{d \\choose \\alpha} = \\frac{d!}{\\alpha_1! \\alpha_2! \\cdots \\alpha_n!}.
```

For instance, create a polynomial with the basis ``[xy^2, xy]`` creates the polynomial
``\\sqrt{3} a xy^2 + \\sqrt{2} b xy`` where `a` and `b` are new JuMP decision variables.
Constraining the polynomial ``axy^2 + bxy`` to be zero with the scaled monomial basis constrains
`a/√3` and `b/√2` to be zero.

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, **2012**.
"""
struct ScaledMonomialBasis{MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
    monomials::MV
end
ScaledMonomialBasis(monomials) = ScaledMonomialBasis(monovec(monomials))

MultivariatePolynomials.polynomialtype(mb::ScaledMonomialBasis{MT}, T::Type) where MT = MultivariatePolynomials.polynomialtype(MT, promote_type(T, Float64))
scaling(m::MultivariatePolynomials.AbstractMonomial) = √(factorial(degree(m)) / prod(factorial, exponents(m)))
MultivariatePolynomials.polynomial(f::Function, mb::ScaledMonomialBasis) = polynomial(i -> scaling(mb.monomials[i]) * f(i), mb.monomials)
unscale_coef(t::MultivariatePolynomials.AbstractTerm) = coefficient(t) / scaling(monomial(t))
function MultivariatePolynomials.coefficients(p, ::Type{<:ScaledMonomialBasis})
    return unscale_coef.(terms(p))
end
