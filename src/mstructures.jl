struct MStruct{B<:AbstractMonomialIndexed,V,E,I,BT<:SA.AbstractBasis{Polynomial{B,V,E},I}} <: SA.MultiplicativeStructure{Polynomial{B,V,E},I}
    basis::BT
end

function (m::MStruct)(a::Polynomial, b::Polynomial, ::Type{U}) where {U}
    return m(m[a], m[b], U)
end

function (m::MStruct{B,V,E})(a::E, b::E, ::Type{Polynomial{B,V,E}}) where {B,V,E}
    return SA.map_keys(Base.Fix1(getindex, m), m(a, b, E))
end
