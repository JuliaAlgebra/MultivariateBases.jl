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
MultivariateBases.Monomial
MultivariateBases.ScaledMonomial
```

## Orthogonal basis

```@docs
MultivariateBases.AbstractMultipleOrthogonal
univariate_orthogonal_basis
reccurence_first_coef
reccurence_second_coef
reccurence_third_coef
reccurence_deno_coef
MultivariateBases.ProbabilistsHermite
MultivariateBases.PhysicistsHermite
MultivariateBases.Laguerre
MultivariateBases.AbstractGegenbauer
MultivariateBases.Legendre
MultivariateBases.ChebyshevFirstKind
MultivariateBases.ChebyshevSecondKind
```
