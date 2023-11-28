"""
    struct OrthonormalCoefficientsBasis{PT<:MP.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis
        polynomials::PV
    end

Polynomial basis with the polynomials of the vector `polynomials` that are
orthonormal with respect to the inner produce derived from the inner product
of their coefficients.
For instance, `FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x])` is the Chebyshev
polynomial basis for cubic polynomials in the variable `x`.
"""
struct OrthonormalCoefficientsBasis{
    PT<:MP.AbstractPolynomialLike,
    PV<:AbstractVector{PT},
} <: AbstractPolynomialVectorBasis{PT,PV}
    polynomials::PV
end

function LinearAlgebra.dot(
    p::MP.AbstractPolynomialLike{S},
    q::MP.AbstractPolynomialLike{T},
    ::Type{<:OrthonormalCoefficientsBasis},
) where {S,T}
    s = zero(MA.promote_operation(*, S, T))
    terms_p = MP.terms(p)
    terms_q = MP.terms(q)
    tsp = iterate(terms_p)
    tsq = iterate(terms_q)
    while !isnothing(tsp) && !isnothing(tsq)
        tp, sp = tsp
        tq, sq = tsq
        cmp = MP.compare(MP.monomial(tp), MP.monomial(tq))
        if iszero(cmp)
            s += conj(MP.coefficient(tp)) * MP.coefficient(tq)
            tsp = iterate(terms_p, sp)
            tsq = iterate(terms_q, sq)
        elseif cmp < 0
            tsp = iterate(terms_p, sp)
        else
            tsq = iterate(terms_q, sq)
        end
    end
    return s
end

function MP.coefficients(p, basis::OrthonormalCoefficientsBasis)
    return [LinearAlgebra.dot(q, p, typeof(basis)) for q in generators(basis)]
end
