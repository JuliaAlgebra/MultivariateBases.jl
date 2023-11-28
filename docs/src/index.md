```@meta
CurrentModule = MultivariateBases
DocTestSetup = quote
    using MultivariateBases
end
```

# MultivariateBases

[MultivariateBases.jl](https://github.com/JuliaAlgebra/MultivariateBases.jl) is a standardized API for multivariate polynomial bases
based on the [MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) API.

```@docs
AbstractPolynomialBasis
maxdegree_basis
basis_covering_monomials
FixedPolynomialBasis
OrthonormalCoefficientsBasis
```

## Monomial basis

```@docs
MonomialBasis
ScaledMonomialBasis
```

## Orthogonal basis

```@docs
AbstractMultipleOrthogonalBasis
univariate_orthogonal_basis
reccurence_first_coef
reccurence_second_coef
reccurence_third_coef
reccurence_deno_coef
ProbabilistsHermiteBasis
PhysicistsHermiteBasis
LaguerreBasis
AbstractGegenbauerBasis
LegendreBasis
ChebyshevBasisFirstKind
ChebyshevBasisSecondKind
```
