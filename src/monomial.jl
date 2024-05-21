abstract type AbstractMonomialIndexed end

struct FullBasis{B<:AbstractMonomialIndexed,M<:MP.AbstractMonomial} <:
       SA.ImplicitBasis{Polynomial{B,M},M} end

function Base.getindex(::FullBasis{B,M}, mono::M) where {B,M}
    return Polynomial{B}(mono)
end

function Base.getindex(::FullBasis{B,M}, p::Polynomial{B,M}) where {B,M}
    return p.monomial
end

SA.mstructure(::FullBasis{B}) where {B} = Mul{B}()

MP.monomial_type(::Type{<:FullBasis{B,M}}) where {B,M} = M
function MP.polynomial_type(basis::FullBasis{B,M}, ::Type{T}) where {B,M,T}
    return MP.polynomial_type(typeof(basis), T)
end

# TODO Move it to SA as :
# struct SubBasis{T,I,B<:SA.ImplicitBasis{T,I},V<:AbstractVector{I}} <: SA.ExplicitBasis{T,Int}
#     implicit::B
#     indices::V
# end
struct SubBasis{B<:AbstractMonomialIndexed,M,V<:AbstractVector{M}} <:
       SA.ExplicitBasis{Polynomial{B,M},Int}
    monomials::V
end

# Overload some of the `AbstractVector` interface for convenience
Base.isempty(basis::SubBasis) = isempty(basis.monomials)
Base.eachindex(basis::SubBasis) = eachindex(basis.monomials)
_iterate(::SubBasis, ::Nothing) = nothing
function _iterate(basis::SubBasis{B}, elem_state) where {B}
    return parent(basis)[elem_state[1]], elem_state[2]
end
Base.iterate(basis::SubBasis) = _iterate(basis, iterate(basis.monomials))
Base.iterate(basis::SubBasis, s) = _iterate(basis, iterate(basis.monomials, s))
Base.length(basis::SubBasis) = length(basis.monomials)
Base.firstindex(basis::SubBasis) = firstindex(basis.monomials)
Base.lastindex(basis::SubBasis) = lastindex(basis.monomials)

MP.nvariables(basis::SubBasis) = MP.nvariables(basis.monomials)
MP.variables(basis::SubBasis) = MP.variables(basis.monomials)

Base.parent(::SubBasis{B,M}) where {B,M} = FullBasis{B,M}()

function Base.getindex(basis::SubBasis, index::Int)
    return parent(basis)[basis.monomials[index]]
end

function Base.getindex(basis::SubBasis{B,M}, value::Polynomial{B,M}) where {B,M}
    j = parent(basis)[value]
    i = searchsortedfirst(basis.indices, j)
    if i in eachindex(basis.indices) && basis.indices[i] == j
        return i
    end
    throw(BoundsError(basis, value))
end

MP.monomial_type(::Type{<:SubBasis{B,M}}) where {B,M} = M

# The `i`th index of output is the index of occurence of `x[i]` in `y`,
# or `0` if it does not occur.
function multi_findsorted(x, y)
    I = zeros(Int, length(x))
    j = 1
    for i in eachindex(x)
        while j ≤ length(y) && x[i] > y[j]
            j += 1
        end
        if j ≤ length(y) && x[i] == y[j]
            I[i] = j
        end
    end
    return I
end

function merge_bases(basis1::MB, basis2::MB) where {MB<:SubBasis}
    monos = MP.merge_monomial_vectors([basis1.monomials, basis2.monomials])
    I1 = multi_findsorted(monos, basis1.monomials)
    I2 = multi_findsorted(monos, basis2.monomials)
    return MB(monos), I1, I2
end

# Unsafe because we don't check that `monomials` is sorted and without duplicates
function unsafe_basis(
    ::Type{B},
    monomials::AbstractVector{M},
) where {B<:AbstractMonomialIndexed,M<:MP.AbstractMonomial}
    return SubBasis{B,M,typeof(monomials)}(monomials)
end

function Base.getindex(::FullBasis{B,M}, monomials::AbstractVector) where {B,M}
    return unsafe_basis(B, MP.monomial_vector(monomials)::AbstractVector{M})
end

function SubBasis{B}(
    monomials::AbstractVector,
) where {B<:AbstractMonomialIndexed}
    return unsafe_basis(
        B,
        MP.monomial_vector(monomials)::AbstractVector{<:MP.AbstractMonomial},
    )
end

function Base.copy(basis::SubBasis)
    return typeof(basis)(copy(basis.monomials))
end

function empty_basis(
    ::Type{<:SubBasis{B,M}},
) where {B<:AbstractMonomialIndexed,M}
    return unsafe_basis(B, MP.empty_monomial_vector(M))
end

function maxdegree_basis(
    ::Type{B},
    variables,
    maxdegree::Int,
) where {B<:AbstractMonomialIndexed}
    return unsafe_basis(B, MP.monomials(variables, 0:maxdegree))
end

function MP.polynomial(coefs::Vector, basis::SubBasis)
    return MP.polynomial(Base.Fix1(getindex, coefs), basis)
end

function MP.polynomial(f::Function, basis::SubBasis)
    if isempty(basis)
        return zero(MP.polynomial_type(basis))
    else
        return MP.polynomial(
            mapreduce(
                ip -> f(ip[1]) * MP.polynomial(ip[2]),
                MA.add!!,
                enumerate(basis),
            ),
        )
    end
end

_convert(::Type{P}, p) where {P} = convert(P, p)
_convert(::Type{P}, ::MA.Zero) where {P} = zero(P)

function MP.polynomial(Q::AbstractMatrix, basis::SubBasis, ::Type{T}) where {T}
    n = length(basis)
    @assert size(Q) == (n, n)
    PT = MP.polynomial_type(eltype(basis), T)
    return _convert(
        PT,
        mapreduce(
            row -> begin
                adjoint(MP.polynomial(basis[row])) * mapreduce(
                    col -> Q[row, col] * MP.polynomial(basis[col]),
                    MA.add!!,
                    1:n;
                    init = MA.Zero(),
                )
            end,
            MA.add!!,
            1:n;
            init = MA.Zero(),
        ),
    )::PT
end

function _show(io::IO, mime::MIME, basis::SubBasis{B}) where {B}
    print(io, "SubBasis{$(nameof(B))}")
    print(io, "([")
    first = true
    # TODO use Base.show_vector here, maybe by wrapping the `generator` vector
    #      into something that spits objects wrapped with the `mime` type
    for mono in basis.monomials
        if !first
            print(io, ", ")
        end
        first = false
        show(io, mime, mono)
    end
    return print(io, "])")
end

function Base.show(io::IO, mime::MIME"text/plain", basis::SubBasis)
    return _show(io, mime, basis)
end

function Base.show(io::IO, mime::MIME"text/print", basis::SubBasis)
    return _show(io, mime, basis)
end

function Base.print(io::IO, basis::SubBasis)
    return show(io, MIME"text/print"(), basis)
end

function Base.show(io::IO, basis::SubBasis)
    return show(io, MIME"text/plain"(), basis)
end

abstract type AbstractMonomial <: AbstractMonomialIndexed end

function basis_covering_monomials(
    ::Type{B},
    monos::AbstractVector,
) where {B<:AbstractMonomial}
    return SubBasis{B}(monos)
end

"""
    struct Monomial <: AbstractMonomialIndexed end

Monomial basis with the monomials of the vector `monomials`.
For instance, `SubBasis{Monomial}([1, x, y, x^2, x*y, y^2])` is the monomial basis
for the subspace of quadratic polynomials in the variables `x`, `y`.

This basis is orthogonal under a scalar product defined with the complex Gaussian measure as density.
Once normalized so as to be orthonormal with this scalar product,
one get ths [`ScaledMonomial`](@ref).
"""
struct Monomial <: AbstractMonomial end

MP.polynomial(p::Polynomial{Monomial}) = MP.polynomial(p.monomial)

function MP.polynomial_type(
    ::Union{SubBasis{B,M},Type{<:SubBasis{B,M}}},
    ::Type{T},
) where {B,M,T}
    return MP.polynomial_type(FullBasis{B,M}, T)
end

function MP.polynomial_type(::Type{FullBasis{B,M}}, ::Type{T}) where {B,M,T}
    return MP.polynomial_type(Polynomial{B,M}, T)
end

function MP.polynomial_type(
    ::Type{Polynomial{Monomial,M}},
    ::Type{T},
) where {M,T}
    return MP.polynomial_type(M, T)
end

function MP.polynomial(f::Function, mb::SubBasis{Monomial})
    return MP.polynomial(f, mb.monomials)
end

function MP.polynomial(Q::AbstractMatrix, mb::SubBasis{Monomial}, T::Type)
    return MP.polynomial(Q, mb.monomials, T)
end

function MP.coefficients(p, basis::SubBasis{Monomial})
    return MP.coefficients(p, basis.monomials)
end

function MP.coefficients(p, ::FullBasis{Monomial})
    return MP.coefficients(p)
end

# Overload some of the `MP` interface for convenience
MP.mindegree(basis::SubBasis{Monomial}) = MP.mindegree(basis.monomials)
MP.maxdegree(basis::SubBasis) = MP.maxdegree(basis.monomials)
MP.extdegree(basis::SubBasis{Monomial}) = MP.extdegree(basis.monomials)
function MP.mindegree(basis::SubBasis{Monomial}, v)
    return MP.mindegree(basis.monomials, v)
end
function MP.maxdegree(basis::SubBasis, v)
    return MP.maxdegree(basis.monomials, v)
end
function MP.extdegree(basis::SubBasis{Monomial}, v)
    return MP.extdegree(basis.monomials, v)
end
