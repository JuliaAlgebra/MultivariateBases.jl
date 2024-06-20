"""
    struct FixedBasis{T} <: SA.ExplicitBasis{T,I}
        elements::Vector{T}
    end

Polynomial basis with the polynomials of the vector `polynomials`.
For instance, `FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x])` is the Chebyshev
polynomial basis for cubic polynomials in the variable `x`.
"""
struct FixedPolynomialBasis{
    PT<:MP.AbstractPolynomialLike,
    PV<:AbstractVector{PT},
} <: AbstractPolynomialVectorBasis{PT,PV}
    polynomials::PV
end
