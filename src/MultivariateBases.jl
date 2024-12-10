module MultivariateBases

import MutableArithmetics as MA
import StarAlgebras as SA
import MultivariatePolynomials as MP

export FullBasis, SubBasis
export maxdegree_basis, explicit_basis_covering, empty_basis, monomial_index
include("interface.jl")

export AbstractMonomialIndexed, Monomial, ScaledMonomial
include("polynomial.jl")
MP.monomial_type(::Type{<:SA.AlgebraElement{A}}) where {A} = MP.monomial_type(A)
function MP.polynomial_type(::Type{<:SA.AlgebraElement{A,T}}) where {A,T}
    return MP.polynomial_type(A, T)
end
struct Algebra{BT,B,M} <:
       SA.AbstractStarAlgebra{Polynomial{B,M},Polynomial{B,M}}
    basis::BT
end
MP.monomial_type(::Type{<:Algebra{B}}) where {B} = MP.monomial_type(B)
function MP.polynomial_type(::Type{<:Algebra{B}}, ::Type{T}) where {B,T}
    return MP.polynomial_type(B, T)
end
function MA.promote_operation(
    ::typeof(SA.basis),
    ::Type{<:Algebra{B}},
) where {B}
    return B
end
SA.basis(a::Algebra) = a.basis

#Base.:(==)(::Algebra{BT1,B1,M}, ::Algebra{BT2,B2,M}) where {BT1,B1,BT2,B2,M} = true
function Base.:(==)(a::Algebra, b::Algebra)
    # `===` is a shortcut for speedup
    return a.basis === b.basis || a.basis == b.basis
end

function Base.show(io::IO, ::Algebra{BT,B}) where {BT,B}
    ioc = IOContext(io, :limit => true, :compact => true)
    return print(ioc, "Polynomial algebra of $B basis")
end

include("monomial.jl")
include("scaled.jl")

export AbstractMultipleOrthogonal,
    ProbabilistsHermite, PhysicistsHermite, Laguerre
export AbstractGegenbauer,
    Legendre, Chebyshev, ChebyshevFirstKind, ChebyshevSecondKind, Trigonometric
export algebra_element,
    sparse_coefficients,
    univariate_orthogonal_basis,
    reccurence_first_coef,
    reccurence_second_coef,
    reccurence_third_coef,
    reccurence_deno_coef

import LinearAlgebra
include("orthogonal.jl")
include("hermite.jl")
include("laguerre.jl")
include("legendre.jl")
include("chebyshev.jl")
include("trigonometric.jl")
include("lagrange.jl")
include("quotient.jl")

function algebra(
    basis::Union{QuotientBasis{Polynomial{B,M}},FullBasis{B,M},SubBasis{B,M}},
) where {B,M}
    return Algebra{typeof(basis),B,M}(basis)
end

function MA.promote_operation(
    ::typeof(algebra),
    BT::Type{
        <:Union{QuotientBasis{Polynomial{B,M}},FullBasis{B,M},SubBasis{B,M}},
    },
) where {B,M}
    return Algebra{BT,B,M}
end

include("arithmetic.jl")
include("fixed.jl")

end # module
