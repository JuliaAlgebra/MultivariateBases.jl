const _APL = MP.AbstractPolynomialLike
# We don't define it for all `AlgebraElement` as this would be type piracy
const _AE = SA.AlgebraElement{<:Algebra}

Base.:(+)(p::_APL, q::_AE) = +(p, MP.polynomial(q))
Base.:(+)(p::_AE, q::_APL) = +(MP.polynomial(p), q)
Base.:(-)(p::_APL, q::_AE) = -(p, MP.polynomial(q))
Base.:(-)(p::_AE, q::_APL) = -(MP.polynomial(p), q)

Base.:(+)(p, q::_AE) = +(constant_algebra_element(typeof(SA.basis(q)), p), q)
Base.:(+)(p::_AE, q) = +(MP.polynomial(p), constant_algebra_element(typeof(SA.basis(p)), q))
Base.:(-)(p, q::_AE) = -(constant_algebra_element(typeof(SA.basis(q)), p), MP.polynomial(q))
Base.:(-)(p::_AE, q) = -(MP.polynomial(p), constant_algebra_element(typeof(SA.basis(p)), q))

function Base.:(+)(p::_AE, q::_AE)
    return MA.operate_to!(SA._preallocate_output(+, p, q), +, p, q)
end

function Base.:(-)(p::_AE, q::_AE)
    return MA.operate_to!(SA._preallocate_output(-, p, q), -, p, q)
end

Base.:(*)(p::Union{_APL,_AE}, q::Polynomial) = *(p, algebra_element(q))
function Base.:(*)(α, p::Polynomial{B}) where {B}
    return _algebra_element(MP.term(α, p.monomial), B)
end