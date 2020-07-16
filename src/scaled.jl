"""
    struct ScaledMonomialBasis{MT<:MP.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
        elements::MV
    end

*Scaled monomial basis* (see [Section 3.1.5, BPT12]) with the monomials of the vector `elements`.
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

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, **2012**.
"""
struct ScaledMonomialBasis{MT<:MP.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractMonomialBasis{MT, MV}
    elements::MV
end
ScaledMonomialBasis(monomials) = ScaledMonomialBasis(monovec(monomials))

MP.polynomialtype(::ScaledMonomialBasis{MT}, T::Type) where MT = MP.polynomialtype(MT, promote_type(T, Float64))
MP.polynomial(f::Function, basis::ScaledMonomialBasis) = MP.polynomial(i -> scaling(basis.elements[i]) * f(i), basis.elements)

function Base.promote_rule(::Type{ScaledMonomialBasis{MT, MV}}, ::Type{MonomialBasis{MT, MV}}) where {MT, MV}
    return MonomialBasis{MT, MV}
end

function change_basis(Q::AbstractMatrix, basis::ScaledMonomialBasis{MT, MV}, B::Type{MonomialBasis{MT, MV}}) where {MT, MV}
    n = length(basis)
    scalings = map(scaling, basis.elements)
    scaled_Q = [Q[i, j] * scalings[i] * scalings[j] for i in 1:n, j in 1:n]
    return scaled_Q, MonomialBasis(basis.elements)
end

function MP.polynomial(Q::AbstractMatrix, basis::ScaledMonomialBasis{MT, MV}, T::Type) where {MT, MV}
    return MP.polynomial(change_basis(Q, basis, MonomialBasis{MT, MV})..., T)
end

scaling(m::MP.AbstractMonomial) = √(factorial(MP.degree(m)) / prod(factorial, MP.exponents(m)))
unscale_coef(t::MP.AbstractTerm) = coefficient(t) / scaling(monomial(t))
function MP.coefficients(p, ::Type{<:ScaledMonomialBasis})
    return unscale_coef.(MP.terms(p))
end
