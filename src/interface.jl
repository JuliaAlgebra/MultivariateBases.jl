"""
    abstract type AbstractPolynomialBasis end

Polynomial basis of a subspace of the polynomials [Section~3.1.5, BPT12].

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, **2012**.
"""
abstract type AbstractPolynomialBasis end

generators(basis::AbstractPolynomialBasis) = basis.polynomials

function Base.copy(basis::AbstractPolynomialBasis)
    return typeof(basis)(copy(generators(basis)))
end
function Base.getindex(
    basis::AbstractPolynomialBasis,
    I::AbstractVector{<:Integer},
)
    return typeof(basis)(generators(basis)[I])
end

# Overload some of the `AbstractVector` interface for convenience
Base.isempty(basis::AbstractPolynomialBasis) = isempty(generators(basis))
Base.eachindex(basis::AbstractPolynomialBasis) = eachindex(generators(basis))
Base.iterate(basis::AbstractPolynomialBasis) = iterate(generators(basis))
Base.iterate(basis::AbstractPolynomialBasis, s) = iterate(generators(basis), s)
Base.length(basis::AbstractPolynomialBasis) = length(generators(basis))
Base.firstindex(basis::AbstractPolynomialBasis) = firstindex(generators(basis))
Base.lastindex(basis::AbstractPolynomialBasis) = lastindex(generators(basis))
Base.getindex(basis::AbstractPolynomialBasis, i::Int) = generators(basis)[i]

# Overload some of the `MP` interface for convenience
MP.mindegree(basis::AbstractPolynomialBasis) = MP.mindegree(generators(basis))
MP.maxdegree(basis::AbstractPolynomialBasis) = MP.maxdegree(generators(basis))
MP.extdegree(basis::AbstractPolynomialBasis) = MP.extdegree(generators(basis))
function MP.mindegree(basis::AbstractPolynomialBasis, v)
    return MP.mindegree(generators(basis), v)
end
function MP.maxdegree(basis::AbstractPolynomialBasis, v)
    return MP.maxdegree(generators(basis), v)
end
function MP.extdegree(basis::AbstractPolynomialBasis, v)
    return MP.extdegree(generators(basis), v)
end
MP.nvariables(basis::AbstractPolynomialBasis) = MP.nvariables(generators(basis))
MP.variables(basis::AbstractPolynomialBasis) = MP.variables(generators(basis))

function MP.polynomial(coefs::Vector, basis::AbstractPolynomialBasis)
    return MP.polynomial(i -> coefs[i], basis)
end

"""
    maxdegree_basis(B::Type{<:AbstractPolynomialBasis}, variables, maxdegree::Int)

Return the basis of type `B` generating all polynomials of degree up to
`maxdegree` with variables `variables`.
"""
function maxdegree_basis end

"""
    basis_covering_monomials(B::Type{<:AbstractPolynomialBasis}, monos::AbstractVector{<:AbstractMonomial})

Return the minimal basis of type `B` that can generate all polynomials of the
monomial basis generated by `monos`.

## Examples

For example, to generate all the polynomials with nonzero coefficients for the
monomials `x^4` and `x^2`, we need three polynomials as otherwise, we generate
polynomials with nonzero constant term.
```jldoctest
julia> using DynamicPolynomials

julia> @polyvar x
(x,)

julia> basis_covering_monomials(ChebyshevBasis, [x^2, x^4])
ChebyshevBasisFirstKind([1.0, -1.0 + 2.0x², 1.0 - 8.0x² + 8.0x⁴])
```
"""
function basis_covering_monomials end

function _show(io::IO, mime::MIME, basis::AbstractPolynomialBasis)
    T = typeof(basis)
    print(io, nameof(T))
    print(io, "([")
    first = true
    # TODO use Base.show_vector here, maybe by wrapping the `generator` vector
    #      into something that spits objects wrapped with the `mime` type
    for g in generators(basis)
        if !first
            print(io, ", ")
        end
        first = false
        show(io, mime, g)
    end
    return print(io, "])")
end

function Base.show(
    io::IO,
    mime::MIME"text/plain",
    basis::AbstractPolynomialBasis,
)
    return _show(io, mime, basis)
end
function Base.show(
    io::IO,
    mime::MIME"text/print",
    basis::AbstractPolynomialBasis,
)
    return _show(io, mime, basis)
end

function Base.print(io::IO, basis::AbstractPolynomialBasis)
    return show(io, MIME"text/print"(), basis)
end
function Base.show(io::IO, basis::AbstractPolynomialBasis)
    return show(io, MIME"text/plain"(), basis)
end
