"""
    struct OrthonormalCoefficientsBasis{T,M<:AbstractMatrix{T},B} <: SA.ExplicitBasis
        matrix::M
        basis::B
    end

Polynomial basis with the polynomials `algebra_element(matrix[i, :], basis)` for
each row `i` of `matrix`. The basis is orthonormal with respect to the inner
product derived from the inner product of their coefficients.
For instance,
```jldoctest
julia> OrthonormalCoefficientsBasis(
    [1 0 0 0
     0 1 0 0
    -1 0 2 0
      -3   4],
    SubBasis{Monomial}([1, x, x^2, x^3]),
)
```
is the Chebyshev polynomial basis for cubic polynomials in the variable `x`.
"""
struct OrthonormalCoefficientsBasis{T,B,M,V} <: SA.ExplicitBasis{
    SA.AlgebraElement{Algebra{SubBasis{B,M,V},B,M},T,Vector{T}},
    Int,
}
    matrix::Matrix{T}
    basis::SubBasis{B,M,V}
end

Base.length(basis::OrthonormalCoefficientsBasis) = size(basis.matrix, 1)

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
