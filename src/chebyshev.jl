abstract type AbstractChebyshevBasis{P} <: AbstractGegenbauerBasis{P} end

function MP.polynomial_type(::Type{<:AbstractChebyshevBasis}, V::Type)
    return MP.polynomial_type(V, Float64)
end

reccurence_first_coef(::Type{<:AbstractChebyshevBasis}, degree) = 2
reccurence_third_coef(::Type{<:AbstractChebyshevBasis}, degree) = -1
reccurence_deno_coef(::Type{<:AbstractChebyshevBasis}, degree) = 1

"""
    struct ChebyshevBasisFirstKind{P} <: AbstractChebyshevBasis{P}
        polynomials::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\frac{1}{\\sqrt{1 - x^2}}`` over the interval ``[-1, 1]``.
"""
struct ChebyshevBasisFirstKind{P} <: AbstractChebyshevBasis{P}
    polynomials::Vector{P}
end

const ChebyshevBasis{P} = ChebyshevBasisFirstKind{P}

function degree_one_univariate_polynomial(
    ::Type{<:ChebyshevBasisFirstKind},
    variable::MP.AbstractVariable,
)
    MA.@rewrite(variable + 0)
end

function scalar_product_function(::Type{<:ChebyshevBasisFirstKind}, i::Int)
    if i == 0
        return π
    elseif isodd(i)
        return 0
    else
        n = div(i, 2)
        return (π / 2^i) * prod(n+1:i) / factorial(n)
    end
end

"""
    struct ChebyshevBasisSecondKind{P} <: AbstractChebyshevBasis{P}
        polynomials::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\sqrt{1 - x^2}`` over the interval ``[-1, 1]``.
"""
struct ChebyshevBasisSecondKind{P} <: AbstractChebyshevBasis{P}
    polynomials::Vector{P}
end

function degree_one_univariate_polynomial(
    ::Type{<:ChebyshevBasisSecondKind},
    variable::MP.AbstractVariable,
)
    MA.@rewrite(2variable + 0)
end

function scalar_product_function(::Type{<:ChebyshevBasisSecondKind}, i::Int)
    if i == 0
        return π / 2
    elseif isodd(i)
        return 0
    else
        n = div(i, 2)
        return π / (2^(i + 1)) * prod(n+2:i) / factorial(n)
    end
end
