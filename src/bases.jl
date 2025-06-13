struct Variables{B,V}
    variables::V
end

function (v::Variables{B})(exponents) where {B}
    return Polynomial{B}(MP.monomial(v.variables, exponents))
end

const FullBasis{B,M} = SA.MappedBasis{Polynomial{B,M}}
const SubBasis{B,M} = SA.SubBasis{Polynomial{B,M}}

function FullBasis{B,M}(vars) where {B,M}
    O = typeof(MP.ordering(vars))
    return SA.MappedBasis{Polynomial{B,M}}(
        MP.ExponentsIterator{O}(vars),
        Variables{B,typeof(vars)}(vars),
        MP.exponents,
    )
end

MP.monomial_type(::Type{<:FullBasis{B,M}}) where {B,M} = M
function MP.polynomial_type(basis::FullBasis{B,M}, ::Type{T}) where {B,M,T}
    return MP.polynomial_type(typeof(basis), T)
end

MP.nvariables(v::Variables) = length(v.variables)
MP.nvariables(basis::FullBasis) = MP.nvariables(basis.map)
MP.nvariables(basis::SubBasis) = MP.nvariables(parent(basis))
MP.variables(v::Variables) = v.variables
MP.variables(basis::FullBasis) = MP.variables(basis.map)
MP.variables(basis::SubBasis) = MP.variables(parent(basis))

function monomial_index(basis::SubBasis, mono::MP.AbstractMonomial)
    @assert MP.variables(basis) == MP.variables(mono)
    return get(basis, MP.exponents(mono), nothing)
end

function explicit_basis_covering(::FullBasis{B}, target::SubBasis{B}) where {B}
    return SubBasis{B}(target.monomials)
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

SubBasis{B}(monos::Tuple) where {B} = SubBasis{B}([monos...])

function Base.copy(basis::SubBasis)
    return typeof(basis)(copy(basis.monomials))
end

function Base.:(==)(a::SubBasis{B}, b::SubBasis{B}) where {B}
    return a.monomials == b.monomials
end

function algebra_type(::Type{BT}) where {B,M,BT<:MonomialIndexedBasis{B,M}}
    return Algebra{BT,B,M}
end

implicit_basis(::SubBasis{B,M}) where {B,M} = FullBasis{B,M}()
implicit_basis(basis::FullBasis) = basis

function implicit(a::SA.AlgebraElement)
    basis = implicit_basis(SA.basis(a))
    return algebra_element(SA.coeffs(a, basis), basis)
end

function MA.promote_operation(
    ::typeof(implicit),
    ::Type{E},
) where {AG,T,E<:SA.AlgebraElement{AG,T}}
    BT = MA.promote_operation(implicit_basis, MA.promote_operation(SA.basis, E))
    A = MA.promote_operation(algebra, BT)
    M = MP.monomial_type(BT)
    return SA.AlgebraElement{A,T,SA.SparseCoefficients{M,T,Vector{M},Vector{T}}}
end

function MA.promote_operation(
    ::typeof(implicit_basis),
    ::Type{<:Union{FullBasis{B,M},SubBasis{B,M}}},
) where {B,M}
    return FullBasis{B,M}
end

function _explicit_basis(coeffs, ::FullBasis{B}) where {B}
    return SubBasis{B}(SA.keys(coeffs))
end

_explicit_basis(_, basis::SubBasis) = basis

function explicit_basis(p::MP.AbstractPolynomialLike)
    return SubBasis{Monomial}(MP.monomials(p))
end

function explicit_basis(a::SA.AlgebraElement)
    return _explicit_basis(SA.coeffs(a), SA.basis(a))
end

function explicit_basis_type(::Type{FullBasis{B,M}}) where {B,M}
    return SubBasis{B,M,MP.monomial_vector_type(M)}
end

function empty_basis(
    ::Type{<:SubBasis{B,M}},
) where {B<:AbstractMonomialIndexed,M}
    return unsafe_basis(B, MP.empty_monomial_vector(M))
end

function maxdegree_basis(
    ::FullBasis{B},
    variables,
    maxdegree::Int,
) where {B<:AbstractMonomialIndexed}
    return unsafe_basis(B, MP.monomials(variables, 0:maxdegree))
end

MP.variables(c::SA.AbstractCoefficients) = MP.variables(SA.keys(c))

_lazy_collect(v::AbstractVector) = collect(v)
_lazy_collect(v::Vector) = v

function sparse_coefficients(p::MP.AbstractPolynomial)
    return SA.SparseCoefficients(
        _lazy_collect(MP.monomials(p)),
        _lazy_collect(MP.coefficients(p)),
    )
end

function sparse_coefficients(t::MP.AbstractTermLike)
    return SA.SparseCoefficients((MP.monomial(t),), (MP.coefficient(t),))
end

function MA.promote_operation(
    ::typeof(sparse_coefficients),
    ::Type{P},
) where {P<:MP.AbstractPolynomialLike}
    M = MP.monomial_type(P)
    T = MP.coefficient_type(P)
    return SA.SparseCoefficients{M,T,Vector{M},Vector{T}}
end

function algebra_element(p::MP.AbstractPolynomialLike)
    return algebra_element(
        sparse_coefficients(p),
        FullBasis{Monomial,MP.monomial_type(p)}(),
    )
end

function algebra_element(f::Function, basis::SubBasis)
    return algebra_element(map(f, eachindex(basis)), basis)
end

_one_if_type(α) = α
_one_if_type(::Type{T}) where {T} = one(T)

function constant_algebra_element_type(
    ::Type{BT},
    ::Type{T},
) where {B,M,BT<:FullBasis{B,M},T}
    A = MA.promote_operation(algebra, BT)
    return SA.AlgebraElement{A,T,SA.SparseCoefficients{M,T,Vector{M},Vector{T}}}
end

function constant_algebra_element(::Type{FullBasis{B,M}}, α) where {B,M}
    return algebra_element(
        sparse_coefficients(
            MP.polynomial(MP.term(_one_if_type(α), MP.constant_monomial(M))),
        ),
        FullBasis{B,M}(),
    )
end

function constant_algebra_element_type(
    ::Type{B},
    ::Type{T},
) where {B<:SubBasis,T}
    A = MA.promote_operation(algebra, B)
    return SA.AlgebraElement{A,T,Vector{T}}
end

function constant_algebra_element(::Type{<:SubBasis{B,M}}, α) where {B,M}
    return algebra_element(
        [_one_if_type(α)],
        SubBasis{B}([MP.constant_monomial(M)]),
    )
end