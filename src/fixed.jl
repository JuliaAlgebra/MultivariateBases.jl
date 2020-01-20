"""
    struct FixedPolynomialBasis{PT<:MP.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis
        polynomials::PV
    end

Polynomial basis with the polynomials of the vector `polynomials`.
For instance, `FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x])` is the Chebyshev
polynomial basis for cubic polynomials in the variable `x`.
"""
struct FixedPolynomialBasis{PT<:MP.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis
    polynomials::PV
end

Base.length(basis::FixedPolynomialBasis) = length(basis.polynomials)
empty_basis(::Type{FixedPolynomialBasis{PT, Vector{PT}}}) where PT = FixedPolynomialBasis(PT[])
function MP.polynomialtype(mb::FixedPolynomialBasis{PT}, T::Type) where PT
    C = MP.coefficienttype(PT)
    U = typeof(zero(C) * zero(T) + zero(C) * zero(T))
    MP.polynomialtype(PT, U)
end
function MP.polynomial(f::Function, fpb::FixedPolynomialBasis)
    sum(ip -> f(ip[1]) * ip[2], enumerate(fpb.polynomials))
end
