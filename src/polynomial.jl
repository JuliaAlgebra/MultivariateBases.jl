# Methods of `Polynomial` for the bases defined in this package. The ones for
# `Monomial` are defined in MultivariatePolynomials.
const _Basis = Union{ScaledMonomial,AbstractMultipleOrthogonal}

function MP.monomial(p::Polynomial{<:_Basis})
    return MP.monomial(MP.variables(p), MP.exponents(p))
end

function Base.hash(p::Polynomial{B}, u::UInt) where {B<:_Basis}
    return hash(
        Tuple(
            (v, e) for
            (v, e) in zip(MP.variables(p), MP.exponents(p)) if !iszero(e)
        ),
        hash(B, u),
    )
end

MP.isconstant(p::Polynomial{<:_Basis}) = all(iszero, MP.exponents(p))

function MP._show(io::IO, mime::MIME, p::Polynomial{B}) where {B<:_Basis}
    print(io, B, "(")
    MP._show(io, mime, MP.monomial(p))
    print(io, ")")
    return
end

function Base.zero(::Type{Polynomial{B,V,E}}) where {B<:_Basis,V,E}
    sc = SA.SparseCoefficients(E[], Int[])
    return SA.algebra_element(sc, algebra(FullBasis{B}(V())))
end

function Base.zero(p::Polynomial{B}) where {B<:_Basis}
    sc = SA.SparseCoefficients(typeof(MP.exponents(p))[], Int[])
    return SA.algebra_element(sc, algebra(FullBasis{B}(p)))
end

function Base.adjoint(p::Polynomial{B}) where {B<:_Basis}
    mono = adjoint(MP.monomial(p))
    return Polynomial(Variables{B}(MP.variables(mono)), MP.exponents(mono))
end

MP.left_constant_mult(α, p::Polynomial{ScaledMonomial}) = SA.Term(α, p)
function MP.right_constant_mult(p::Polynomial{ScaledMonomial}, α)
    return MP.left_constant_mult(α, p)
end
