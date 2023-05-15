"""
    struct ScaledMonomialBasis{MT<:MP.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
        monomials::MV
    end

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
struct ScaledMonomialBasis{MT<:MP.AbstractMonomial,MV<:AbstractVector{MT}} <:
       AbstractMonomialBasis{MT,MV}
    monomials::MV
end
function ScaledMonomialBasis(monomials)
    return ScaledMonomialBasis(MP.monomial_vector(monomials))
end

function Base.getindex(basis::ScaledMonomialBasis, i::Int)
    mono = basis.monomials[i]
    return scaling(mono) * mono
end

function MP.polynomial_type(::ScaledMonomialBasis{MT}, T::Type) where {MT}
    return MP.polynomial_type(MT, promote_type(T, Float64))
end
function MP.polynomial(f::Function, basis::ScaledMonomialBasis)
    return MP.polynomial(
        i -> scaling(basis.monomials[i]) * f(i),
        basis.monomials,
    )
end

function Base.promote_rule(
    ::Type{ScaledMonomialBasis{MT,MV}},
    ::Type{MonomialBasis{MT,MV}},
) where {MT,MV}
    return MonomialBasis{MT,MV}
end

function change_basis(
    Q::AbstractMatrix,
    basis::ScaledMonomialBasis{MT,MV},
    B::Type{MonomialBasis{MT,MV}},
) where {MT,MV}
    n = length(basis)
    scalings = map(scaling, basis.monomials)
    scaled_Q = [Q[i, j] * scalings[i] * scalings[j] for i in 1:n, j in 1:n]
    return scaled_Q, MonomialBasis(basis.monomials)
end

function MP.polynomial(
    Q::AbstractMatrix,
    basis::ScaledMonomialBasis{MT,MV},
    T::Type,
) where {MT,MV}
    return MP.polynomial(change_basis(Q, basis, MonomialBasis{MT,MV})..., T)
end

function scaling(m::MP.AbstractMonomial)
    return √(factorial(MP.degree(m)) / prod(factorial, MP.exponents(m)))
end
unscale_coef(t::MP.AbstractTerm) = MP.coefficient(t) / scaling(MP.monomial(t))
function MP.coefficients(p, ::Type{<:ScaledMonomialBasis})
    return unscale_coef.(MP.terms(p))
end
