module MultivariateBases

import MutableArithmetics as MA

import MultivariatePolynomials as MP

export AbstractPolynomialBasis
export maxdegree_basis, basis_covering_monomials, empty_basis
include("interface.jl")

export AbstractMonomialBasis, MonomialBasis, ScaledMonomialBasis
include("monomial.jl")
include("scaled.jl")

export FixedPolynomialBasis,
    AbstractMultipleOrthogonalBasis,
    ProbabilistsHermiteBasis,
    PhysicistsHermiteBasis,
    LaguerreBasis
export AbstractGegenbauerBasis,
    LegendreBasis,
    ChebyshevBasis,
    ChebyshevBasisFirstKind,
    ChebyshevBasisSecondKind
export generators,
    univariate_orthogonal_basis,
    reccurence_first_coef,
    reccurence_second_coef,
    reccurence_third_coef,
    reccurence_deno_coef
include("fixed.jl")
include("orthogonal.jl")
include("hermite.jl")
include("laguerre.jl")
include("legendre.jl")
include("chebyshev.jl")

end # module
