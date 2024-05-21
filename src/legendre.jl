"""
    struct AbstractGegenbauer <: AbstractMultipleOrthogonal end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = (1 - x^2)^{\\alpha - 1/2}`` over the interval ``[-1, 1]``.
"""
abstract type AbstractGegenbauer <: AbstractMultipleOrthogonal end

even_odd_separated(::Type{<:AbstractGegenbauer}) = true
reccurence_second_coef(::Type{<:AbstractGegenbauer}, degree) = 0

"""
    struct Legendre <: AbstractGegenbauer end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = 1`` over the interval ``[-1, 1]``.
"""
struct Legendre <: AbstractGegenbauer end

function MP.polynomial_type(
    ::Type{Polynomial{Legendre,M}},
    ::Type{T},
) where {M,T}
    return MP.polynomial_type(M, float(T))
end

reccurence_first_coef(::Type{Legendre}, degree) = (2degree - 1)
reccurence_third_coef(::Type{Legendre}, degree) = -(degree - 1)
reccurence_deno_coef(::Type{Legendre}, degree) = degree

function degree_one_univariate_polynomial(
    ::Type{Legendre},
    variable::MP.AbstractVariable,
)
    MA.@rewrite(variable + 0)
end

function _scalar_product_function(::Type{Legendre}, i::Int)
    if isodd(i)
        return 0
    else
        return 2 / (i + 1)
    end
end
