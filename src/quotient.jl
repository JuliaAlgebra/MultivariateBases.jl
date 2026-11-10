struct QuotientBasis{T,I,B<:SA.AbstractBasis{T,I},D} <: SA.ExplicitBasis{T,I}
    basis::B
    divisor::D
end

implicit_basis(basis::QuotientBasis) = implicit_basis(basis.basis)
MP.variables(basis::QuotientBasis) = MP.variables(basis.basis)

Base.iterate(basis::QuotientBasis) = iterate(basis.basis)
Base.iterate(basis::QuotientBasis, s) = iterate(basis.basis, s)
Base.length(basis::QuotientBasis) = length(basis.basis)

function SA.algebra_element(coeffs, basis::QuotientBasis)
    return SA.AlgebraElement(coeffs, algebra(basis))
end

function MP.coefficients(p, basis::QuotientBasis)
    return MP.coefficients(rem(p, basis.divisor), basis.basis)
end

function SA.coeffs(coeffs, b::FullBasis{Monomial}, basis::QuotientBasis)
    return MP.coefficients(
        MP.polynomial(values(coeffs), keys_as_monomials(keys(coeffs), b)),
        basis,
    )
end

function SA.coeffs(coeffs, sub::SubBasis{Monomial}, basis::QuotientBasis)
    return MP.coefficients(MP.polynomial(coeffs, keys_as_monomials(sub)), basis)
end

function SA.adjoint_coeffs(
    coeffs,
    src::SubBasis{Monomial},
    dest::QuotientBasis{<:BasisElement{Monomial}},
)
    return map(src) do poly
        return sum(MP.terms(rem(MP.polynomial(poly), dest.divisor))) do t
            return MP.coefficient(t) *
                   coeffs[SA.key_index(dest.basis, MP.exponents(t))]
        end
    end
end

"""
    struct LinearQuotientBasis{T,I,B<:SA.AbstractBasis{T,I},M<:AbstractMatrix} <:
        SA.ExplicitBasis{T,I}
        basis::B
        matrix::M
    end

Basis of the quotient of the span of `basis` by the kernel of
`transpose(matrix)`. The coefficients of a polynomial `p` in this basis are
`transpose(matrix) * a` where `a` is the vector of coefficients of `p` in
`basis`. In particular, `p` is zero in the quotient if and only if its vector
of coefficients is zero.

For instance, if `matrix` has orthonormal columns spanning the orthogonal
complement of a subspace `K` of the span of `basis`, a polynomial is zero in
this basis if and only if it belongs to `K`. Compared to `QuotientBasis`,
this does not require a `rem` function but the polynomial must belong to the
span of `basis`.
"""
struct LinearQuotientBasis{T,I,B<:SA.AbstractBasis{T,I},M<:AbstractMatrix} <:
       SA.ExplicitBasis{T,I}
    basis::B
    matrix::M
    function LinearQuotientBasis(
        basis::SA.AbstractBasis{T,I},
        matrix::AbstractMatrix,
    ) where {T,I}
        if size(matrix, 1) != length(basis)
            throw(
                DimensionMismatch(
                    "The number of rows `$(size(matrix, 1))` of the matrix does not match the length `$(length(basis))` of the basis",
                ),
            )
        end
        return new{T,I,typeof(basis),typeof(matrix)}(basis, matrix)
    end
end

implicit_basis(basis::LinearQuotientBasis) = implicit_basis(basis.basis)
MP.variables(basis::LinearQuotientBasis) = MP.variables(basis.basis)

Base.length(basis::LinearQuotientBasis) = size(basis.matrix, 2)

function Base.show(io::IO, basis::LinearQuotientBasis)
    print(
        io,
        "Quotient of dimension $(length(basis)) of the span of ",
        basis.basis,
    )
    return
end

function SA.algebra_element(coeffs, basis::LinearQuotientBasis)
    return SA.AlgebraElement(coeffs, algebra(basis))
end

function _key_index(basis::SubBasis{Monomial}, mono::MP.AbstractMonomialLike)
    return SA.key_index(
        basis,
        MP.exponents(MP.monomial(mono), MP.variables(basis)),
    )
end

function _linear_quotient_coeffs(p, basis::LinearQuotientBasis)
    sub = basis.basis
    for mono in MP.monomials(p)
        if isnothing(_key_index(sub, mono))
            throw(
                ArgumentError(
                    "The monomial `$mono` is not in the basis `$sub` of which `$basis` is a quotient.",
                ),
            )
        end
    end
    return _transpose_mul(basis.matrix, MP.coefficients(p, sub))
end

# Computes `transpose(W) * a` skipping the zero entries of `W`. This matters
# when the elements of `a` are, e.g., affine functions of decision variables
# as the product with zero would still create a term with zero coefficient.
function _transpose_mul(W::AbstractMatrix, a::AbstractVector)
    U = MA.promote_operation(*, eltype(W), eltype(a))
    return map(axes(W, 2)) do k
        acc = zero(U)
        for i in axes(W, 1)
            if !iszero(W[i, k])
                acc = MA.operate!!(MA.add_mul, acc, W[i, k], a[i])
            end
        end
        return acc
    end
end

function MP.coefficients(p, basis::LinearQuotientBasis)
    return _linear_quotient_coeffs(p, basis)
end

function SA.coeffs(coeffs, b::FullBasis{Monomial}, basis::LinearQuotientBasis)
    return _linear_quotient_coeffs(
        MP.polynomial(values(coeffs), keys_as_monomials(keys(coeffs), b)),
        basis,
    )
end

function SA.coeffs(coeffs, sub::SubBasis{Monomial}, basis::LinearQuotientBasis)
    return _linear_quotient_coeffs(
        MP.polynomial(coeffs, keys_as_monomials(sub)),
        basis,
    )
end

function SA.adjoint_coeffs(
    coeffs,
    src::SubBasis{Monomial},
    dest::LinearQuotientBasis{<:BasisElement{Monomial}},
)
    a = dest.matrix * coeffs
    return map(keys_as_monomials(src)) do mono
        i = _key_index(dest.basis, mono)
        return isnothing(i) ? zero(eltype(a)) : a[i]
    end
end
