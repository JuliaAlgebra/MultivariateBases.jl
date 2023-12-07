abstract type AbstractHermite <: AbstractMultipleOrthogonal end

function MP.polynomial_type(::Type{Polynomial{B,M}}, ::Type{T}) where {B<:AbstractHermite,M,T}
    return MP.polynomial_type(M, float(T))
end

even_odd_separated(::Type{<:AbstractHermite}) = true

reccurence_second_coef(::Type{<:AbstractHermite}, degree) = 0
reccurence_deno_coef(::Type{<:AbstractHermite}, degree) = 1

"""
    struct ProbabilistsHermiteBasis{P} <: AbstractHermiteBasis{P}
        polynomials::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\exp(-x^2/2)`` over the interval ``[-\\infty, \\infty]``.
"""
struct ProbabilistsHermite <: AbstractHermite end
reccurence_first_coef(::Type{ProbabilistsHermite}, degree) = 1
function reccurence_third_coef(::Type{ProbabilistsHermite}, degree)
    return -(degree - 1)
end
function degree_one_univariate_polynomial(
    ::Type{ProbabilistsHermite},
    variable::MP.AbstractVariable,
)
    MA.@rewrite(1variable)
end

function _scalar_product_function(::Type{ProbabilistsHermite}, i::Int)
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
    struct PhysicistsHermite{P} <: AbstractHermite{P}
        polynomials::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\exp(-x^2)`` over the interval ``[-\\infty, \\infty]``.
"""
struct PhysicistsHermite <: AbstractHermite end
reccurence_first_coef(::Type{PhysicistsHermite}, degree) = 2
reccurence_third_coef(::Type{PhysicistsHermite}, degree) = -2(degree - 1)
function degree_one_univariate_polynomial(
    ::Type{PhysicistsHermite},
    variable::MP.AbstractVariable,
)
    MA.@rewrite(2variable)
end

function _scalar_product_function(::Type{PhysicistsHermite}, i::Int)
    if i == 0
        return √(π)
    elseif isodd(i)
        return 0
    else
        n = div(i, 2)
        return (√(π) / (2^i)) * prod(n+1:i)
    end
end
