# Multivariate Bases

[![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] [![Build Status][build-img]][build-url] [![Codecov branch][codecov-img]][codecov-url]

This package provides a standardized API for multivariate polynomial bases
based on the [MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) API.

It defines the following basis:
* `FixedBasis`: A polynomial basis described by a list of polynomials.
* Monomial bases: `MonomialBasis` and `ScaledMonomialBasis`.
* Orthogonal bases:
  - Hermite bases: `ProbabilistsHermiteBasis` and `PhysicistsHermiteBasis`.
  - `LaguerreBasis`.
  - Gegenbauer bases:
    * `LegendreBasis`.
    * Chebyshev bases: `ChebyshevBasisFirstKind` and `ChebyshevBasisSecondKind`

See the documentation for more details.

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://JuliaAlgebra.github.io/MultivariateBases.jl/stable
[docs-latest-url]: https://JuliaAlgebra.github.io/MultivariateBases.jl/dev

[build-img]: https://github.com/JuliaAlgebra/MultivariateBases.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/JuliaAlgebra/MultivariateBases.jl/actions?query=workflow%3ACI
[codecov-img]: http://codecov.io/github/JuliaAlgebra/MultivariateBases.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaAlgebra/MultivariateBases.jl?branch=master
