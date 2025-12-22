struct MStruct{
    B<:AbstractMonomialIndexed,
    V,
    E,
    I,
    BT<:SA.AbstractBasis{Polynomial{B,V,E},I},
} <: SA.MultiplicativeStructure{Polynomial{B,V,E},I}
    basis::BT
end

function Base.:(==)(a::MStruct, b::MStruct)
    return a.basis == b.basis
end

function MA.promote_operation(
    ::typeof(SA.basis),
    ::Type{MStruct{B,V,E,I,BT}},
) where {B,V,E,I,BT}
    return BT
end

function (m::MStruct)(a::Polynomial, b::Polynomial, ::Type{U}) where {U}
    return m(m[a], m[b], U)
end

function (m::MStruct{B,V,E})(
    a::E,
    b::E,
    ::Type{Polynomial{B,V,E}},
) where {B,V,E}
    return SA.map_keys(Base.Fix1(getindex, m), m(a, b, E))
end

SA.promote_with_map(::MStruct, b, m) = MStruct(b), m

function SA.promote_basis_with_maps(a::SA.AbstractBasis, b::MStruct)
    _a, _b = SA.promote_basis_with_maps(a, SA.basis(b))
    return _a, SA.maybe_promote(b, _b...)
end

SA.promote_object(v::Variables, m::MStruct, map) = SA.promote_object(v, SA.basis(m), map)
SA.promote_object(v::Variables, m::SA.SubBasis, map) = SA.promote_object(v, parent(m), map)
SA.promote_object(::Variables, m::FullBasis, _) = m.map
