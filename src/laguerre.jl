"""
    struct LaguerreBasis <: AbstractMultipleOrthogonal end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\exp(-x)`` over the interval ``[0, \\infty]``.
"""
struct Laguerre <: AbstractMultipleOrthogonal end

# TODO implement multiplication with https://www.jstor.org/stable/2002985

_promote_coef(::Type{T}, ::Type{<:Laguerre}) where {T} = _float(T)

even_odd_separated(::Type{Laguerre}) = false

reccurence_first_coef(::Type{Laguerre}, degree) = -1
reccurence_second_coef(::Type{Laguerre}, degree) = (2degree - 1)
reccurence_third_coef(::Type{Laguerre}, degree) = -(degree - 1)
reccurence_deno_coef(::Type{Laguerre}, degree) = degree

function degree_one_univariate_polynomial(::Type{Laguerre}, variable)
    MA.@rewrite(1 - variable)
end

function _scalar_product_function(::Type{Laguerre}, i::Int)
    return factorial(i)
end
