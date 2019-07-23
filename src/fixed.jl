"""
    struct FixedPolynomialBasis{PT<:MultivariatePolynomials.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis
        polynomials::PV
    end

Polynomial basis with the polynomials of the vector `polynomials`.
For instance, `FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x])` is the Chebyshev
polynomial basis for cubic polynomials in the variable `x`.
"""
struct FixedPolynomialBasis{PT<:MultivariatePolynomials.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis
    polynomials::PV
end

function MultivariatePolynomials.polynomialtype(mb::FixedPolynomialBasis{PT}, T::Type) where PT
    C = MultivariatePolynomials.coefficienttype(PT)
    U = typeof(zero(C) * zero(T) + zero(C) * zero(T))
    MultivariatePolynomials.polynomialtype(PT, U)
end
function MultivariatePolynomials.polynomial(f::Function, fpb::FixedPolynomialBasis)
    sum(ip -> f(ip[1]) * ip[2], enumerate(fpb.polynomials))
end
