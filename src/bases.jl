const FullBasis{B,V,E} = SA.MappedBasis{Polynomial{B,V,E}}
const SubBasis{B,V,E} = SA.SubBasis{Polynomial{B,V,E}}
const MonomialIndexedBasis{B,V,E} = Union{SubBasis{B,V,E},FullBasis{B,V,E}}

function FullBasis{B}(vars) where {B}
    O = typeof(MP.ordering(vars))
    v = Variables{B,typeof(vars)}(vars)
    exps = MP.ExponentsIterator{O}(constant_monomial_exponents(v))
    return SA.MappedBasis{Polynomial{B,typeof(vars),eltype(exps)}}(exps, v, MP.exponents)
end

function FullBasis{B}(p::MP.AbstractPolynomialLike) where {B}
    return FullBasis{B}(MP.variables(p))
end

SA.star(b::MonomialIndexedBasis{B,V,E}, exp::E) where {B,V,E} = b[SA.star(b[exp])]

constant_monomial_exponents(b::FullBasis) = constant_monomial_exponents(b.map)

MP.monomial_type(::Type{<:FullBasis{B,V}}) where {B,V} = MP.monomial_type(V)
function MP.polynomial_type(basis::FullBasis, ::Type{T}) where {T}
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

function explicit_basis_covering(full::FullBasis{B}, target::SubBasis{B}) where {B}
    return SA.SubBasis(full, target.keys)
end

MP.monomial_type(::Type{<:SA.SubBasis{T,I,K,B}}) where {T,I,K,B} = MP.monomial_type(B)

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
    monomials::AbstractVector{<:MP.AbstractMonomial},
) where {B<:AbstractMonomialIndexed}
    # TODO We should add a way to directly get the vector of exponents inside `DynamicPolynomials.MonomialVector`.
    # `MP.exponents.(monomials)` should work in the meantime even if it's not the most efficient
    return SA.SubBasis(FullBasis{B}(MP.variables(monomials)), MP.exponents.(monomials))
end

function Base.getindex(::FullBasis{B}, monomials::AbstractVector{M}) where {B,M<:MP.AbstractMonomial}
    return unsafe_basis(B, MP.monomial_vector(monomials)::AbstractVector{M})
end

function SubBasis{B}(
    monomials::AbstractVector{<:MP.AbstractTermLike},
) where {B<:AbstractMonomialIndexed}
    return unsafe_basis(
        B,
        MP.monomial_vector(monomials)::AbstractVector{<:MP.AbstractMonomial},
    )
end

SubBasis{B}(monos::Tuple) where {B} = SubBasis{B}([monos...])

function algebra_type(::Type{BT}) where {B,M,BT<:MonomialIndexedBasis{B,M}}
    return Algebra{BT,B,M}
end

implicit_basis(basis::SubBasis) = parent(basis)
implicit_basis(basis::FullBasis) = basis

function implicit(a::SA.AlgebraElement)
    basis = implicit_basis(SA.basis(a))
    return algebra_element(SA.coeffs(a, basis), basis)
end

function MA.promote_operation(
    ::typeof(implicit),
    ::Type{AE},
) where {AG,T,AE<:SA.AlgebraElement{AG,T}}
    BT = MA.promote_operation(implicit_basis, MA.promote_operation(SA.basis, AE))
    A = MA.promote_operation(algebra, BT)
    E = SA.key_type(BT)
    return SA.AlgebraElement{A,T,SA.SparseCoefficients{E,T,Vector{E},Vector{T},typeof(isless)}}
end

MA.promote_operation(::typeof(implicit_basis), B::Type{<:SA.ImplicitBasis}) = B
MA.promote_operation(::typeof(implicit_basis), ::Type{<:SA.SubBasis{T,I,K,B}}) where {T,I,K,B} = B

function _explicit_basis(coeffs, basis::FullBasis{B}) where {B}
    return SA.SubBasis(basis, _lazy_collect(SA.keys(coeffs)))
end

_explicit_basis(_, basis::SubBasis) = basis

function explicit_basis(p::MP.AbstractPolynomialLike)
    return SubBasis{Monomial}(MP.monomials(p))
end

function explicit_basis(a::SA.AlgebraElement)
    return _explicit_basis(SA.coeffs(a), SA.basis(a))
end

function explicit_basis_type(BT::Type{<:FullBasis{B,V,E}}) where {B,V,E}
    return SA.SubBasis{eltype(BT),Int,SA.key_type(BT),BT,Vector{E}}
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
_lazy_collect(v::Tuple) = collect(v)
_lazy_collect(v::Vector) = v

function sparse_coefficients(p::MP.AbstractPolynomial)
    return SA.SparseCoefficients(
        MP.exponents.(MP.monomials(p)),
        _lazy_collect(MP.coefficients(p)),
    )
end

function sparse_coefficients(t::MP.AbstractTermLike)
    return SA.SparseCoefficients((MP.exponents(t),), (MP.coefficient(t),))
end

function algebra_element(p::MP.AbstractPolynomialLike)
    return algebra_element(
        sparse_coefficients(p),
        FullBasis{Monomial}(p),
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
) where {B,V,E,BT<:FullBasis{B,V,E},T}
    A = MA.promote_operation(algebra, BT)
    return SA.AlgebraElement{A,T,SA.SparseCoefficients{E,T,Tuple{E},Tuple{T},typeof(isless)}}
end

function constant_algebra_element(b::FullBasis, α)
    return algebra_element(
        SA.SparseCoefficients((constant_monomial_exponents(b),), (_one_if_type(α),)),
        b,
    )
end

function constant_algebra_element_type(
    ::Type{B},
    ::Type{T},
) where {B<:SubBasis,T}
    A = MA.promote_operation(algebra, B)
    return SA.AlgebraElement{A,T,Vector{T}}
end

function constant_algebra_element(basis::SubBasis, α)
    return algebra_element(
        [_one_if_type(α)],
        SA.SubBasis(parent(basis), [constant_monomial_exponents(parent(basis))]),
    )
end
