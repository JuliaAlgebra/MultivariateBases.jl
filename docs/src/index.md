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
maxdegree_basis
explicit_basis_covering
```

## Basis elements

```@docs
Polynomial
SemisimpleElement
```

## Monomial basis

```@docs
Monomial
ScaledMonomial
```

## Orthogonal basis

```@docs
AbstractMultipleOrthogonal
univariate_orthogonal_basis
reccurence_first_coef
reccurence_second_coef
reccurence_third_coef
reccurence_deno_coef
ProbabilistsHermite
PhysicistsHermite
Laguerre
AbstractGegenbauer
Legendre
ChebyshevFirstKind
ChebyshevSecondKind
Trigonometric
```

## Additional basis

```@docs
FixedBasis
SemisimpleBasis
```
