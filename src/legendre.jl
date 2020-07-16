"""
    struct AbstractGegenbauerBasis{P} <: AbstractMultipleOrthogonalBasis{P}
        elements::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = (1 - x^2)^{\\alpha - 1/2}`` over the interval ``[-1, 1]``.
"""
abstract type AbstractGegenbauerBasis{P} <: AbstractMultipleOrthogonalBasis{P} end

even_odd_separated(::Type{<:AbstractGegenbauerBasis}) = true
reccurence_second_coef(::Type{<:AbstractGegenbauerBasis}, degree) = 0

"""
    struct LegendreBasis{P} <: AbstractGegenbauerBasis{P}
        elements::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = 1`` over the interval ``[-1, 1]``.
"""
struct LegendreBasis{P} <: AbstractGegenbauerBasis{P}
    elements::Vector{P}
end

polynomial_type(::Type{<:LegendreBasis}, V::Type) = MP.polynomialtype(V, Float64)

reccurence_first_coef(::Type{<:LegendreBasis}, degree) = (2degree - 1)
reccurence_third_coef(::Type{<:LegendreBasis}, degree) = -(degree - 1)
reccurence_deno_coef(::Type{<:LegendreBasis}, degree) = degree

degree_one_univariate_polynomial(::Type{<:LegendreBasis}, variable::MP.AbstractVariable) = MA.@rewrite(variable + 0)

function scalar_product_function(::Type{<:LegendreBasis}, i::Int)
    if isodd(i)
        return 0
    else
        return 2/(i+1)
    end
end
