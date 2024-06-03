module MultivariateBases

import MutableArithmetics as MA
import StarAlgebras as SA
import MultivariatePolynomials as MP

export AbstractPolynomialBasis, FullBasis, SubBasis
export maxdegree_basis, basis_covering_monomials, empty_basis
include("interface.jl")

export AbstractMonomialIndexed, Monomial, ScaledMonomial
include("polynomial.jl")
MP.monomial_type(::Type{<:SA.AlgebraElement{A}}) where {A} = MP.monomial_type(A)
const Algebra{BT,B,M} = SA.StarAlgebra{Polynomial{B,M},Polynomial{B,M},BT}
MP.monomial_type(::Type{<:Algebra{B}}) where {B} = MP.monomial_type(B)
include("monomial.jl")
include("scaled.jl")

export AbstractMultipleOrthogonal,
    ProbabilistsHermite, PhysicistsHermite, Laguerre
export AbstractGegenbauer,
    Legendre, Chebyshev, ChebyshevFirstKind, ChebyshevSecondKind
export generators,
    univariate_orthogonal_basis,
    reccurence_first_coef,
    reccurence_second_coef,
    reccurence_third_coef,
    reccurence_deno_coef
#include("fixed.jl")

import LinearAlgebra
#include("orthonormal.jl")
include("orthogonal.jl")
include("hermite.jl")
include("laguerre.jl")
include("legendre.jl")
include("chebyshev.jl")
include("quotient.jl")

SA.algebra(basis::Union{QuotientBasis,FullBasis,SubBasis}) = SA.StarAlgebra(_object(basis), basis)

end # module
