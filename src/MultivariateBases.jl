module MultivariateBases

import MutableArithmetics
const MA = MutableArithmetics

using MultivariatePolynomials
const MP = MultivariatePolynomials

export AbstractPolynomialBasis
export maxdegree_basis, empty_basis
include("interface.jl")

export AbstractMonomialBasis, MonomialBasis, ScaledMonomialBasis
include("monomial.jl")
include("scaled.jl")

export FixedPolynomialBasis, ChebyshevBasis
include("fixed.jl")
include("chebyshev.jl")

end # module
