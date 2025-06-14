module MultivariateBases

import MutableArithmetics as MA
import StarAlgebras as SA
import MultivariatePolynomials as MP

export FullBasis, SubBasis
export maxdegree_basis, explicit_basis_covering, empty_basis, monomial_index
include("interface.jl")

struct Variables{B,V}
    variables::V
end

Variables{B}(vars) where {B} = Variables{B,typeof(vars)}(vars)

MP.monomial_type(::Type{Variables{B,V}}) where {B,V} = MP.monomial_type(V)

constant_monomial_exponents(v::Variables) = map(_ -> 0, v.variables)

function (v::Variables)(exponents)
    return Polynomial(v, exponents)
end

export AbstractMonomialIndexed, Monomial, ScaledMonomial
include("polynomial.jl")
MP.monomial_type(::Type{<:SA.AlgebraElement{A}}) where {A} = MP.monomial_type(A)
function MP.polynomial_type(::Type{<:SA.AlgebraElement{A,T}}) where {A,T}
    return MP.polynomial_type(A, T)
end
MP.monomial_type(::Type{<:SA.StarAlgebra{O}}) where {O} = MP.monomial_type(O)
function MP.polynomial_type(::Type{A}, ::Type{T}) where {A<:SA.StarAlgebra,T}
    return MP.polynomial_type(MP.monomial_type(A), T)
end

include("bases.jl")
include("mstructures.jl")
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
    return SA.StarAlgebra(Variables{B}(MP.variables(basis)), MStruct(basis))
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
