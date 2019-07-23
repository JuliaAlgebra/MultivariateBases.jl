module MultivariateBases

using MultivariatePolynomials

export AbstractPolynomialBasis
export FixedPolynomialBasis, ScaledMonomialBasis, MonomialBasis

include("interface.jl")

include("monomial.jl")
include("fixed.jl")
include("scaled.jl")

end # module
