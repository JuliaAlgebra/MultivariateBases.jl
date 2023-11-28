abstract type AbstractHermiteBasis{P} <: AbstractMultipleOrthogonalBasis{P} end

function MP.polynomial_type(::Type{<:AbstractHermiteBasis}, V::Type)
    return MP.polynomial_type(V, Int)
end

even_odd_separated(::Type{<:AbstractHermiteBasis}) = true

reccurence_second_coef(::Type{<:AbstractHermiteBasis}, degree) = 0
reccurence_deno_coef(::Type{<:AbstractHermiteBasis}, degree) = 1

"""
    struct ProbabilistsHermiteBasis{P} <: AbstractHermiteBasis{P}
        polynomials::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\exp(-x^2/2)`` over the interval ``[-\\infty, \\infty]``.
"""
struct ProbabilistsHermiteBasis{P} <: AbstractHermiteBasis{P}
    polynomials::Vector{P}
end
reccurence_first_coef(::Type{<:ProbabilistsHermiteBasis}, degree) = 1
function reccurence_third_coef(::Type{<:ProbabilistsHermiteBasis}, degree)
    return -(degree - 1)
end
function degree_one_univariate_polynomial(
    ::Type{<:ProbabilistsHermiteBasis},
    variable::MP.AbstractVariable,
)
    MA.@rewrite(1variable)
end

function scalar_product_function(::Type{<:ProbabilistsHermiteBasis}, i::Int)
    if i == 0
        return √(2 * π)
    elseif isodd(i)
        return 0
    else
        n = div(i, 2)
        return (√(2 * π) / (2^n)) * prod(n+1:2*n)
    end
end

"""
    struct PhysicistsHermiteBasis{P} <: AbstractHermiteBasis{P}
        polynomials::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\exp(-x^2)`` over the interval ``[-\\infty, \\infty]``.
"""
struct PhysicistsHermiteBasis{P} <: AbstractHermiteBasis{P}
    polynomials::Vector{P}
end
reccurence_first_coef(::Type{<:PhysicistsHermiteBasis}, degree) = 2
reccurence_third_coef(::Type{<:PhysicistsHermiteBasis}, degree) = -2(degree - 1)
function degree_one_univariate_polynomial(
    ::Type{<:PhysicistsHermiteBasis},
    variable::MP.AbstractVariable,
)
    MA.@rewrite(2variable)
end

function scalar_product_function(::Type{<:PhysicistsHermiteBasis}, i::Int)
    if i == 0
        return √(π)
    elseif isodd(i)
        return 0
    else
        n = div(i, 2)
        return (√(π) / (2^i)) * prod(n+1:i)
    end
end
