abstract type AbstractHermiteBasis{P} <: AbstractMultipleOrthogonalBasis{P} end

polynomial_type(::Type{<:AbstractHermiteBasis}, V::Type) = MP.polynomialtype(V, Int)

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
reccurence_third_coef(::Type{<:ProbabilistsHermiteBasis}, degree) = -(degree - 1)
degree_one_univariate_polynomial(::Type{<:ProbabilistsHermiteBasis}, variable::MP.AbstractVariable) = MA.@rewrite(1variable)

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
degree_one_univariate_polynomial(::Type{<:PhysicistsHermiteBasis}, variable::MP.AbstractVariable) = MA.@rewrite(2variable)
