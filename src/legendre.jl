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

_promote_coef(::Type{T}, ::Type{<:Legendre}) where {T} = _float(T)

reccurence_first_coef(::Type{Legendre}, degree) = (2degree - 1)
reccurence_third_coef(::Type{Legendre}, degree) = -(degree - 1)
reccurence_deno_coef(::Type{Legendre}, degree) = degree

function degree_one_univariate_polynomial(::Type{Legendre}, variable)
    MA.@rewrite(variable + 0)
end

function _scalar_product_function(::Type{Legendre}, i::Int)
    if isodd(i)
        return 0
    else
        return 2 / (i + 1)
    end
end
