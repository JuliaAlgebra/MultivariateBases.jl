abstract type AbstractMonomial <: AbstractMonomialIndexed end

function explicit_basis_covering(
    ::FullBasis{B},
    target::SubBasis{<:AbstractMonomial},
) where {B<:AbstractMonomial}
    return SubBasis{B}(target.monomials)
end

# To break ambiguity
function explicit_basis_covering(
    ::FullBasis{B},
    target::SubBasis{B},
) where {B<:AbstractMonomial}
    return SubBasis{B}(target.monomials)
end

function Base.adjoint(p::Polynomial{B}) where {B<:AbstractMonomialIndexed}
    mono = adjoint(MP.monomial(p))
    return Polynomial(Variables{B}(MP.variables(mono)), MP.exponents(mono))
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

degree_one_univariate_polynomial(::Type{Monomial}, value) = value

function recurrence_eval(::Type{Monomial}, previous, value, degree)
    return previous[degree] * value
end

function (m::MStruct{Monomial,V,E})(a::E, b::E, ::Type{E}) where {V,E}
    return SA.SparseCoefficients((a .+ b,), (1,))
end

SA.coeffs(p::Polynomial{Monomial}, ::FullBasis{Monomial}) = p.monomial

function MP.polynomial_type(
    ::Union{SubBasis{B,V,E},Type{<:SubBasis{B,V,E}}},
    ::Type{T},
) where {B,V,E,T}
    return _polynomial_type(B, V, T)
end

function MP.polynomial_type(
    ::Type{Polynomial{B,V,E}},
    ::Type{T},
) where {B,V,E,T}
    return _polynomial_type(B, V, T)
end

function MP.polynomial(f::Function, mb::SubBasis{Monomial})
    return MP.polynomial(f, mb.monomials)
end

function MP.polynomial(Q::AbstractMatrix, mb::SubBasis{Monomial}, T::Type)
    return MP.polynomial(Q, mb.monomials, T)
end

function MP.coefficients(
    p::MP.AbstractPolynomialLike,
    basis::SubBasis{Monomial},
)
    return MP.coefficients(p, basis.monomials)
end

function MP.coefficients(p::MP.AbstractPolynomialLike, ::FullBasis{Monomial})
    return p
end

function _assert_constant(α) end

function _assert_constant(
    x::Union{Polynomial,SA.AlgebraElement,MP.AbstractPolynomialLike},
)
    return error("Expected constant element, got type `$(typeof(x))`")
end

#function MA.operate!(::SA.UnsafeAddMul{<:Mul{Monomial}}, p::MP.AbstractPolynomial, args::Vararg{Any,N}) where {N}
#    return MA.operate!(MA.add_mul, p, args...)
#end

#function MA.operate!(
#    ::SA.UnsafeAddMul{Mul{B}},
#    res::SA.AbstractCoefficients,
#    v::SA.AbstractCoefficients,
#    w::SA.AbstractCoefficients,
#) where {B<:MB.AbstractMonomial}
#    for (kv, a) in nonzero_pairs(v)
#        for (kw, b) in nonzero_pairs(w)
#            SA.unsafe_push!(res, kv * kw, a * b)
#            c = ms.structure(kv, kw)
#            for (k, v) in nonzero_pairs(c)
#            end
#        end
#    end
#    return res
#end

#function MA.operate!(op::SA.UnsafeAddMul{Mul{B}}, a::SA.AbstractCoefficients, α, b::SA.AbstractCoefficients) where {B<:MB.AbstractMonomial}
#    for
#    MA.operate!(op, a, α, Polynomial{Monomial}(x.monomial * y.monomial * z.monomial))
#    return a
#end

#function MA.operate!(
#    ::SA.UnsafeAddMul{typeof(*)},
#    a::SA.AlgebraElement{<:SA.StarAlgebra{<:MonomialIndexedBasis,Monomial}},
#    α,
#    x::Polynomial{Monomial},
#)
#    _assert_constant(α)
#    SA.unsafe_push!(a, x, α)
#    return a
#end

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

_promote_coef(::Type{T}, ::Type{Monomial}) where {T} = T

function MP.polynomial_type(::Type{<:FullBasis{B,V}}, ::Type{T}) where {T,B,V}
    return _polynomial_type(B, V, T)
end

function _polynomial_type(::Type{B}, ::Type{V}, ::Type{T}) where {B,V,T}
    return MP.polynomial_type(MP.monomial_type(V), _promote_coef(T, B))
end

_vec(v::Vector) = v
_vec(v::AbstractVector) = collect(v)

# Adapted from SA to incorporate `_promote_coef`
function SA.coeffs(
    cfs,
    source::MonomialIndexedBasis{B1},
    target::MonomialIndexedBasis{B2},
) where {B1,B2}
    source === target && return cfs
    source == target && return cfs
    if B1 === B2 && target isa FullBasis
        # The defaults initialize to zero and then sums which promotes
        # `JuMP.VariableRef` to `JuMP.AffExpr`
        return SA.SparseCoefficients(_vec(source.keys), _vec(cfs))
    else
        res = SA.zero_coeffs(
            _promote_coef(_promote_coef(SA.value_type(cfs), B1), B2),
            target,
        )
        return SA.coeffs!(res, cfs, source, target)
    end
end

# FIXME this assumes that the basis is invariant under adjoint
SA.star(::SubBasis, coeffs) = SA.star.(coeffs)

# TODO use Base.show_vector here, maybe by wrapping the `generator` vector
#      into something that spits objects wrapped with the `mime` type
function _show_vector(io::IO, mime::MIME, v, map = identity)
    print(io, '[')
    first = true
    for el in v
        if !first
            print(io, ", ")
        end
        first = false
        show(io, mime, map(el))
    end
    return print(io, ']')
end

function _show(io::IO, mime::MIME, basis::SubBasis{B}) where {B}
    print(io, "SubBasis{$(nameof(B))}(")
    _show_vector(
        io,
        mime,
        basis.keys,
        Base.Fix1(MP.monomial, MP.variables(basis)),
    )
    print(io, ')')
    return
end

function Base.show(
    io::IO,
    mime::MIME"text/plain",
    basis::Union{Variables,SubBasis},
)
    return _show(io, mime, basis)
end

function Base.show(
    io::IO,
    mime::MIME"text/print",
    basis::Union{Variables,SubBasis},
)
    return _show(io, mime, basis)
end

function Base.print(io::IO, basis::Union{Variables,SubBasis})
    return show(io, MIME"text/print"(), basis)
end

function Base.show(io::IO, basis::Union{Variables,SubBasis})
    return show(io, MIME"text/plain"(), basis)
end
