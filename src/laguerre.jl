"""
    struct LaguerreBasis{P} <: AbstractMultipleOrthogonalBasis{P}
        elements::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\exp(-x)`` over the interval ``[0, \\infty]``.
"""
struct LaguerreBasis{P} <: AbstractMultipleOrthogonalBasis{P}
    elements::Vector{P}
end

polynomial_type(::Type{<:LaguerreBasis}, V::Type) = MP.polynomialtype(V, Float64)

even_odd_separated(::Type{<:LaguerreBasis}) = false

reccurence_first_coef(::Type{<:LaguerreBasis}, degree) = -1
reccurence_second_coef(::Type{<:LaguerreBasis}, degree) = (2degree - 1)
reccurence_third_coef(::Type{<:LaguerreBasis}, degree) = -(degree - 1)
reccurence_deno_coef(::Type{<:LaguerreBasis}, degree) = degree

degree_one_univariate_polynomial(::Type{<:LaguerreBasis}, variable::MP.AbstractVariable) = MA.@rewrite(1 - variable)

function scalar_product_function(::Type{<:LaguerreBasis}, i::Int)
    return factorial(i) 
end

