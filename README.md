# Multivariate Bases

[![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] [![Build Status][build-img]][build-url] [![Codecov branch][codecov-img]][codecov-url]

This package provides a standardized API for multivariate polynomial bases
based on the [MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) API.

This package defines a few univariate bases such as chebyshev polynomials.
These are then extended to multivariate as follows.
Given a monomial `x^i * y^j` the corresponding multivariate chebyshev polynomial is the product of the `i`th univariate chebyshev polynomial in `x` and the `j`th univariate chebyshev polynomial in `y`.
Given this one-to-one correspondence between monomials and multivariate chebyshev polynomials, we represent them directly by these monomials and keep this representation even through addition, multiplication, etc...
```julia
julia> using DynamicPolynomials

julia> @polyvar x y;

julia> using MultivariateBases

julia> basis = FullBasis{Chebyshev,typeof(x*y)}();

julia> basis[x^2 * y^3]
ChebyshevFirstKind(x²y³)

julia> basis[x^2 * y^3] * basis[x * y]
1//4·ChebyshevFirstKind(xy²) + 1//4·ChebyshevFirstKind(xy⁴) + 1//4·ChebyshevFirstKind(x³y²) + 1//4·ChebyshevFirstKind(x³y⁴)
```

The elements obtained by manipulating these polynomials are `StarAlgebras.AlgebraElement`.
The algebra in [StarAlgebras.jl](https://github.com/JuliaAlgebra/StarAlgebras.jl) implements the [MutableArithmetics.jl](https://github.com/jump-dev/MutableArithmetics.jl/) API for efficient manipulation.

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
