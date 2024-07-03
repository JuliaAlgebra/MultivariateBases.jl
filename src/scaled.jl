"""
    struct ScaledMonomial <: AbstractMonomial end

*Scaled monomial basis* (see [Section 3.1.5, BPT12]) with the monomials of the vector `monomials`.
Given a monomial ``x^\\alpha = x_1^{\\alpha_1} \\cdots x_n^{\\alpha_n}`` of degree ``d = \\sum_{i=1}^n \\alpha_i``,
the corresponding polynomial of the basis is
```math
{d \\choose \\alpha}^{\\frac{1}{2}} x^{\\alpha} \\quad \\text{ where } \\quad
{d \\choose \\alpha} = \\frac{d!}{\\alpha_1! \\alpha_2! \\cdots \\alpha_n!}.
```

For instance, create a polynomial with the basis ``[xy^2, xy]`` creates the polynomial
``\\sqrt{3} a xy^2 + \\sqrt{2} b xy`` where `a` and `b` are new JuMP decision variables.
Constraining the polynomial ``axy^2 + bxy`` to be zero with the scaled monomial basis constrains
`a/√3` and `b/√2` to be zero.

This basis is orthonormal under the scalar product:
```math
\\langle f, g \\rangle = \\int_{\\mathcal{C}^n} f(z) \\overline{g(z)} d\\nu_n
```
where ``\\nu_n`` is the Gaussian measure on ``\\mathcal{C}^n`` with the density
``\\pi^{-n} \\exp(-\\lVert z \\rVert^2)``.
See [Section 4; B07] for more details.


[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics (2012).

[B07] Barvinok, Alexander.
*Integration and optimization of multivariate polynomials by restriction onto a random subspace.*
Foundations of Computational Mathematics 7.2 (2007): 229-244.
"""
struct ScaledMonomial <: AbstractMonomial end

function (::Mul{ScaledMonomial})(a::MP.AbstractMonomial, b::MP.AbstractMonomial)
    mono = a * b
    α = prod(
        MP.variables(mono);
        init = inv(binomial(MP.degree(mono), MP.degree(a))),
    ) do v
        return binomial(MP.degree(mono, v), MP.degree(a, v))
    end
    return sparse_coefficients(MP.term(√α, mono))
end

function SA.coeffs(p::Polynomial{ScaledMonomial}, ::FullBasis{Monomial})
    return scaling(p.monomial) * p.monomial
end

_float(::Type{T}) where {T<:Number} = float(T)
# Could be for instance `MathOptInterface.ScalarAffineFunction{Float64}`
# which does not implement `float`
_float(::Type{F}) where {F} = MA.promote_operation(+, F, F)

_promote_coef(::Type{T}, ::Type{ScaledMonomial}) where {T} = _float(T)

function MP.polynomial(f::Function, basis::SubBasis{ScaledMonomial}, ::Type{T}) where {T}
    return MP.polynomial(
        i -> scaling(T, basis.monomials[i]) * f(i),
        basis.monomials,
    )
end

function Base.promote_rule(
    ::Type{SubBasis{ScaledMonomial,M,V}},
    ::Type{SubBasis{Monomial,M,V}},
) where {M,V}
    return SubBasis{Monomial,M,V}
end

function scaling(::Type{T}, m::MP.AbstractMonomial) where {T}
    return √(T(
        factorial(MP.degree(m)) / prod(factorial, MP.exponents(m); init = 1),
    ))
end
unscale_coef(t::MP.AbstractTerm) = MP.coefficient(t) / scaling(MP.coefficient_type(t), MP.monomial(t))
function SA.coeffs(
    t::MP.AbstractTermLike,
    ::FullBasis{ScaledMonomial},
    ::FullBasis{Monomial},
)
    mono = MP.monomial(t)
    return MP.term(mono * MP.coefficient(t), mono)
end
function MP.coefficients(p, ::FullBasis{ScaledMonomial})
    return unscale_coef.(MP.coefficient_type(p), MP.terms(p))
end
function MP.coefficients(p, basis::SubBasis{ScaledMonomial})
    return MP.coefficients(p, basis.monomials) ./ scaling.(MP.coefficient_type(p), MP.monomials(p))
end
function SA.coeffs(
    p::MP.AbstractPolynomialLike,
    ::FullBasis{Monomial},
    basis::FullBasis{ScaledMonomial},
)
    return MP.polynomial(MP.coefficients(p, basis), MP.monomials(p))
end

function SA.coeffs!(
    res,
    cfs,
    source::MonomialIndexedBasis{ScaledMonomial},
    target::MonomialIndexedBasis{Monomial},
)
    MA.operate!(zero, res)
    for (k, v) in SA.nonzero_pairs(cfs)
        mono = source[k].monomial
        SA.unsafe_push!(
            res,
            target[Polynomial{Monomial}(mono)],
            v * scaling(eltype(res), mono),
        )
    end
    MA.operate!(SA.canonical, res)
    return res
end

function SA.coeffs!(
    res,
    cfs,
    source::MonomialIndexedBasis{Monomial},
    target::MonomialIndexedBasis{ScaledMonomial},
)
    MA.operate!(zero, res)
    for (k, v) in SA.nonzero_pairs(cfs)
        mono = source[k].monomial
        SA.unsafe_push!(
            res,
            target[Polynomial{ScaledMonomial}(mono)],
            v / scaling(eltype(T), mono),
        )
    end
    MA.operate!(SA.canonical, res)
    return res
end
