module MultivariateBases

import MutableArithmetics as MA
import StarAlgebras as SA
import MultivariatePolynomials as MP

import MultivariatePolynomials:
    AbstractMonomialIndexed,
    AbstractMonomialBasis,
    Monomial,
    BasisElement,
    Variables,
    FullBasis,
    SubBasis,
    MonomialIndexedBasis,
    MStruct,
    algebra,
    algebra_element,
    sparse_coefficients,
    convert_basis,
    implicit,
    implicit_basis,
    explicit_basis,
    explicit_basis_type,
    explicit_basis_covering,
    empty_basis,
    maxdegree_basis,
    keys_as_monomials,
    exponents_index,
    variable_index,
    constant_algebra_element,
    constant_algebra_element_type,
    algebra_element_type,
    recurrence_eval,
    degree_one_univariate_polynomial,
    _polynomial_type,
    _promote_coef

export FullBasis, SubBasis
export maxdegree_basis, explicit_basis_covering, empty_basis
export keys_as_monomials
export AbstractMonomialIndexed, Monomial, ScaledMonomial
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

function algebra(basis::QuotientBasis{<:BasisElement{B,V}}) where {B,V}
    mstruct = if MP.is_commutative(V)
        MStruct(basis)
    else
        if B != Monomial
            error(
                "Only the `Monomial` basis is supported with noncommutative variables",
            )
        end
        SA.DiracMStructure(basis, *)
    end
    return SA.StarAlgebra(Variables{B}(MP.variables(basis)), mstruct)
end

function MA.promote_operation(
    ::typeof(algebra),
    BT::Type{<:QuotientBasis{BasisElement{B,V,E}}},
) where {B,V,E}
    P = BasisElement{B,V,E}
    MS = if MP.is_commutative(V)
        MStruct{B,V,E,SA.key_type(BT),BT}
    else
        SA.DiracMStructure{P,SA.key_type(BT),BT,typeof(*)}
    end
    return SA.StarAlgebra{Variables{B,V},P,MS}
end

include("arithmetic.jl")
include("fixed.jl")
end # module
