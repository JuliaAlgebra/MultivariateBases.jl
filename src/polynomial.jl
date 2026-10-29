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
MP.monomial_type(::Type{P}) where {P<:Polynomial{<:_Basis}} = P
MP.monomial_type(p::Polynomial{<:_Basis}) = typeof(p)
MP.nvariables(p::Polynomial{<:_Basis}) = MP.nvariables(p.variables)

function MP.algebra_element(p::Polynomial{B}) where {B<:_Basis}
    return algebra_element(
        sparse_coefficients(MP.monomial(p)),
        FullBasis{B}(MP.variables(p)),
    )
end
MP.polynomial(p::Polynomial{<:_Basis}) = algebra_element(p)

function MP._show(io::IO, mime::MIME, p::Polynomial{B}) where {B<:_Basis}
    print(io, B, "(")
    MP._show(io, mime, MP.monomial(p))
    print(io, ")")
    return
end

function Base.show(io::IO, mime::MIME"text/latex", p::Polynomial{<:_Basis})
    print(io, "\$\$ ")
    MP._show(io, mime, p)
    print(io, " \$\$")
    return
end

function Base.show(io::IO, mime::MIME"text/plain", p::Polynomial{<:_Basis})
    return MP._show(io, mime, p)
end

function Base.show(io::IO, mime::MIME"text/print", p::Polynomial{<:_Basis})
    return MP._show(io, mime, p)
end

Base.show(io::IO, p::Polynomial{<:_Basis}) = show(io, MIME"text/plain"(), p)
Base.print(io::IO, p::Polynomial{<:_Basis}) = show(io, MIME"text/print"(), p)

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

MP.term(p::Polynomial{ScaledMonomial}) = SA.Term(true, p)
MP.left_constant_mult(α, p::Polynomial{ScaledMonomial}) = SA.Term(α, p)
function MP.right_constant_mult(p::Polynomial{ScaledMonomial}, α)
    return MP.left_constant_mult(α, p)
end

for op in [:+, :-, :*]
    @eval begin
        function Base.$op(p::Polynomial{<:_Basis}, q::MP.AbstractPolynomial)
            return $op(SA.Term(1, p), q)
        end
        function Base.$op(p::MP.AbstractPolynomial, q::Polynomial{<:_Basis})
            return $op(p, SA.Term(1, q))
        end
    end
end

Base.:*(α::Number, p::Polynomial{<:_Basis}) = SA.algebra_element(SA.Term(α, p))
