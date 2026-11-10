"""
    struct QuotientBasis{T,I,B<:SA.AbstractBasis{T,I},M<:AbstractMatrix} <:
        SA.ExplicitBasis{T,I}
        basis::B
        matrix::M
    end

Basis of the quotient of the span of `basis` by the kernel of `matrix`. The
coefficients of a polynomial `p` in this basis are `matrix * a` where `a` is
the vector of coefficients of `p` in `basis`. In particular, `p` is zero in the
quotient if and only if `matrix * a` is zero. The adjoint map sends the
coefficients `c` in this basis to the coefficients `transpose(matrix) * c` in
`basis`.

To represent the quotient of the span `V` of `basis` modulo an ideal `J`, the
only requirement on `matrix` is that its kernel is `J ∩ V`, i.e., `matrix * a`
is zero if and only if the polynomial with coefficients `a` belongs to `J`.
Recombining the rows of `matrix` with an invertible matrix does not change its
kernel so it represents the same quotient, only its coordinates change. If the
rows of `matrix` are linearly independent, the length of the basis is the
dimension of the quotient. Two typical choices are:

* If a Gröbner basis of `J` is known, `matrix[i, j]` is the coefficient of the
  `i`th standard monomial in the remainder of the `j`th element of `basis`
  modulo the Gröbner basis. For instance, modulo `x^2 + 1`, the basis
  `[1, x, x^2, x^3]` gives `matrix = [1 0 -1 0; 0 1 0 -1]`.
* Otherwise, `J` can be replaced by its truncation `J_d` where `d` is the
  maximal degree of `basis`. `J_d` is the row space of the Macaulay matrix of
  degree `d`, i.e., the span of `m * h` for each generator `h` of `J` and
  monomial `m` such that the degree of `m * h` is at most `d`. The rows of
  `matrix` can be the reduced row echelon form of a basis of the linear
  equations defining `J_d ∩ V`. Then `matrix * a` gives the coefficients of
  the *truncated normal form*: the unique polynomial supported on the pivot
  monomials that is congruent to `p` modulo `J_d`. Note that `J_d ∩ V` may be
  strictly contained in `J ∩ V` if a combination of generators has a degree
  smaller than the generators that are combined (a degree fall).
"""
struct QuotientBasis{T,I,B<:SA.AbstractBasis{T,I},M<:AbstractMatrix} <:
       SA.ExplicitBasis{T,I}
    basis::B
    matrix::M
    function QuotientBasis(
        basis::SA.AbstractBasis{T,I},
        matrix::AbstractMatrix,
    ) where {T,I}
        if size(matrix, 2) != length(basis)
            throw(
                DimensionMismatch(
                    "The number of columns `$(size(matrix, 2))` of the matrix does not match the length `$(length(basis))` of the basis",
                ),
            )
        end
        return new{T,I,typeof(basis),typeof(matrix)}(basis, matrix)
    end
end

implicit_basis(basis::QuotientBasis) = implicit_basis(basis.basis)
MP.variables(basis::QuotientBasis) = MP.variables(basis.basis)

Base.length(basis::QuotientBasis) = size(basis.matrix, 1)

function Base.show(io::IO, basis::QuotientBasis)
    print(
        io,
        "Quotient of dimension $(length(basis)) of the span of ",
        basis.basis,
    )
    return
end

function SA.algebra_element(coeffs, basis::QuotientBasis)
    return SA.AlgebraElement(coeffs, algebra(basis))
end

function _key_index(basis::SubBasis{Monomial}, mono::MP.AbstractMonomialLike)
    return SA.key_index(
        basis,
        MP.exponents(MP.monomial(mono), MP.variables(basis)),
    )
end

function _quotient_coeffs(p, basis::QuotientBasis)
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
    return _mul(basis.matrix, MP.coefficients(p, sub))
end

# Computes `A * a` skipping the zero entries of `A`. This matters when the
# elements of `a` are, e.g., affine functions of decision variables as the
# product with zero would still create a term with zero coefficient.
function _mul(A::AbstractMatrix, a::AbstractVector)
    U = MA.promote_operation(*, eltype(A), eltype(a))
    c = [zero(U) for _ in axes(A, 1)]
    for j in axes(A, 2)
        for i in axes(A, 1)
            if !iszero(A[i, j])
                c[i] = MA.operate!!(MA.add_mul, c[i], A[i, j], a[j])
            end
        end
    end
    return c
end

function _mul(A::SparseArrays.AbstractSparseMatrixCSC, a::AbstractVector)
    U = MA.promote_operation(*, eltype(A), eltype(a))
    c = [zero(U) for _ in axes(A, 1)]
    rows = SparseArrays.rowvals(A)
    vals = SparseArrays.nonzeros(A)
    for j in axes(A, 2)
        for k in SparseArrays.nzrange(A, j)
            if !iszero(vals[k])
                i = rows[k]
                c[i] = MA.operate!!(MA.add_mul, c[i], vals[k], a[j])
            end
        end
    end
    return c
end

function MP.coefficients(p, basis::QuotientBasis)
    return _quotient_coeffs(p, basis)
end

function SA.coeffs(coeffs, b::FullBasis{Monomial}, basis::QuotientBasis)
    return _quotient_coeffs(
        MP.polynomial(values(coeffs), keys_as_monomials(keys(coeffs), b)),
        basis,
    )
end

function SA.coeffs(coeffs, sub::SubBasis{Monomial}, basis::QuotientBasis)
    return _quotient_coeffs(
        MP.polynomial(coeffs, keys_as_monomials(sub)),
        basis,
    )
end

function SA.adjoint_coeffs(
    coeffs,
    src::SubBasis{Monomial},
    dest::QuotientBasis{<:BasisElement{Monomial}},
)
    a = transpose(dest.matrix) * coeffs
    return map(keys_as_monomials(src)) do mono
        i = _key_index(dest.basis, mono)
        return isnothing(i) ? zero(eltype(a)) : a[i]
    end
end
