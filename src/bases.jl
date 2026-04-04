const FullBasis{B,V,E} = SA.MappedBasis{Polynomial{B,V,E}}
const SubBasis{B,V,E} = SA.SubBasis{Polynomial{B,V,E}}
const MonomialIndexedBasis{B,V,E} = Union{SubBasis{B,V,E},FullBasis{B,V,E}}

exponents_index(::FullBasis, exp) = exp
exponents_index(b::SubBasis, exp) = SA.key_index(b, exp)

function MP.ordering(::Type{MonomialIndexedBasis{B,V,E}}) where {B,V,E}
    return MP.ordering(E)
end
MP.ordering(b::MonomialIndexedBasis) = MP.ordering(typeof(b))

# Base.Fix2 will implement `==` by checking `===` on the list of variables
struct _ExponentsWithVariables{V}
    variables::V
end

(e::_ExponentsWithVariables)(mono) = MP.exponents(mono, e.variables)

function Base.:(==)(a::_ExponentsWithVariables, b::_ExponentsWithVariables)
    unique(sort(a.variables)) == unique(sort(b.variables))
end

function FullBasis(vars::Variables{B,V}) where {B,V}
    O = MP.ordering(vars.variables)
    exps = MP.ExponentsIterator{O}(constant_monomial_exponents(vars))
    inverse_map = if MP.is_commutative(V)
        MP.exponents
    else
        _ExponentsWithVariables(vars.variables)
    end
    return SA.MappedBasis{Polynomial{B,V,eltype(exps)}}(
        exps,
        vars,
        inverse_map,
    )
end

function FullBasis{B}(vars) where {B}
    return FullBasis(Variables{B,typeof(vars)}(vars))
end

function FullBasis{B}(p::MP.AbstractPolynomialLike) where {B}
    return FullBasis{B}(MP.variables(p))
end

function SA.star(b::MonomialIndexedBasis{B,V,E}, exp::E) where {B,V,E}
    return b[SA.star(b[exp])]
end

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

function explicit_basis_covering(
    full::FullBasis{B},
    target::SubBasis{B},
) where {B}
    return SA.SubBasis(full, target.keys)
end

function MP.monomial_type(::Type{<:SA.SubBasis{T,I,K,B}}) where {T,I,K,B}
    return MP.monomial_type(B)
end

# FIXME workaround for TP, we should redirect to typeof
MP.monomial_type(b::FullBasis) = MP.monomial_type(b.map)
MP.monomial_type(b::SubBasis) = MP.monomial_type(parent(b))

# FIXME type piracy
SA.comparable(::MP.ExponentsIterator{M}) where {M} = M()

# Unsafe because we don't check that `monomials` is sorted and without duplicates
function unsafe_basis(
    ::Type{B},
    monomials::AbstractVector{<:MP.AbstractMonomial},
) where {B<:AbstractMonomialIndexed}
    # TODO We should add a way to directly get the vector of exponents inside `DynamicPolynomials.MonomialVector`.
    # `MP.exponents.(monomials)` should work in the meantime even if it's not the most efficient
    return SA.SubBasis(
        FullBasis{B}(MP.variables(monomials)),
        MP.exponents.(monomials),
    )
end

function Base.getindex(
    ::FullBasis{B},
    monomials::AbstractVector{M},
) where {B,M<:MP.AbstractMonomial}
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

function _coeffs_type(::Type{C}, ::Type{B}) where {C,B<:FullBasis}
    T = eltype(C) # Even works for `NTuple`!
    E = SA.key_type(B)
    return SA.SparseCoefficients{E,T,_similar_type(C, E),C,typeof(isless)}
end

function _coeffs_type(::Type{C}, ::Type{B}) where {C,B<:SubBasis}
    return C
end

function algebra_element_type(
    ::Type{C}, # type of coefficient vector
    ::Type{B}, # basis type
) where {C,B}
    A = MA.promote_operation(algebra, B)
    return SA.AlgebraElement{
        eltype(C), # Even works for `NTuple`!
        A,
        _coeffs_type(C, B),
    }
end

function MA.promote_operation(
    ::typeof(implicit),
    ::Type{AE},
) where {AG,T,AE<:SA.AlgebraElement{T,AG}}
    BT =
        MA.promote_operation(implicit_basis, MA.promote_operation(SA.basis, AE))
    return algebra_element_type(Vector{T}, BT)
end

MA.promote_operation(::typeof(implicit_basis), B::Type{<:SA.ImplicitBasis}) = B
function MA.promote_operation(
    ::typeof(implicit_basis),
    ::Type{<:SA.SubBasis{T,I,K,B}},
) where {T,I,K,B}
    return B
end

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
    basis::FullBasis{B},
    maxdegree::Int,
) where {B<:AbstractMonomialIndexed}
    return unsafe_basis(B, MP.monomials(MP.variables(basis), 0:maxdegree))
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
    return algebra_element(sparse_coefficients(p), FullBasis{Monomial}(p))
end

function algebra_element(f::Function, basis::SubBasis)
    return algebra_element(map(f, eachindex(basis)), basis)
end

# Another option is to use `fieldcount`, see
# https://discourse.julialang.org/t/get-tuple-length-from-type/32483/16?u=blegat
_similar_type(::Type{<:NTuple{N,Any}}, ::Type{T}) where {N,T} = NTuple{N,T}
function _similar_type(::Type{V}, ::Type{T}) where {V<:AbstractVector,T}
    return SA.similar_type(V, T)
end

function full_basis_type(
    ::Type{B},
    ::Type{P},
) where {B,P<:MP.AbstractPolynomialLike}
    V = MA.promote_operation(MP.variables, P)
    E = _similar_type(V, Int)
    O = MP.ordering(P)
    return SA.MappedBasis{
        Polynomial{B,V,E},
        E,
        MP.ExponentsIterator{O,Nothing,E},
        Variables{B,V},
        typeof(MP.exponents),
    }
end

function MA.promote_operation(
    ::typeof(algebra_element),
    P::Type{<:MP.AbstractPolynomialLike{T}},
) where {T}
    return algebra_element_type(Vector{T}, full_basis_type(Monomial, P))
end

_one_if_type(α) = α
_one_if_type(::Type{T}) where {T} = one(T)

function constant_algebra_element_type(
    ::Type{BT},
    ::Type{T},
) where {BT<:FullBasis,T}
    return algebra_element_type(Tuple{T}, BT)
end

function constant_algebra_element(b::FullBasis, α)
    return algebra_element(
        SA.SparseCoefficients(
            (constant_monomial_exponents(b),),
            (_one_if_type(α),),
        ),
        b,
    )
end

function constant_algebra_element_type(
    ::Type{B},
    ::Type{T},
) where {B<:SubBasis,T}
    A = MA.promote_operation(algebra, B)
    return SA.AlgebraElement{T,A,Vector{T}}
end

function constant_algebra_element(basis::SubBasis, α)
    return algebra_element(
        [_one_if_type(α)],
        SA.SubBasis(
            parent(basis),
            [constant_monomial_exponents(parent(basis))],
        ),
    )
end

function SA.promote_with_map(::Variables{B}, vars, map) where {B}
    return Variables{B}(vars), map
end

SA.promote_with_map(::FullBasis, vars, m) = FullBasis(vars), m

function MP.promote_variables_with_maps(a::Variables, b::Variables)
    _a, _b = MP.promote_variables_with_maps(a.variables, b.variables)
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end

function SA.promote_bases_with_maps(a::FullBasis, b::FullBasis)
    _a, _b = MP.promote_variables_with_maps(a.map, b.map)
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end
