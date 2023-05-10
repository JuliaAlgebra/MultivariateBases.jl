abstract type AbstractChebyshevBasis{P} <: AbstractGegenbauerBasis{P} end

MP.polynomial_type(::Type{<:AbstractChebyshevBasis}, V::Type) = MP.polynomial_type(V, Float64)

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

degree_one_univariate_polynomial(::Type{<:ChebyshevBasisFirstKind}, variable::MP.AbstractVariable) = MA.@rewrite(variable + 0)

"""
    struct ChebyshevBasisSecondKind{P} <: AbstractChebyshevBasis{P}
        polynomials::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\sqrt{1 - x^2}`` over the interval ``[-1, 1]``.
"""
struct ChebyshevBasisSecondKind{P} <: AbstractChebyshevBasis{P}
    polynomials::Vector{P}
end

degree_one_univariate_polynomial(::Type{<:ChebyshevBasisSecondKind}, variable::MP.AbstractVariable) = MA.@rewrite(2variable + 0)
