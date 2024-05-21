# TODO Add to MultivariatePolynomials
Base.keytype(p::MP.AbstractPolynomial) = MP.monomial_type(p)
Base.valtype(p::MP.AbstractPolynomial) = MP.coefficient_type(p)
#Base.keys(p::MP.AbstractPolynomial) = MP.monomials(p)
SA.nonzero_pairs(p::MP.AbstractPolynomial) = MP.terms(p)
function Base.similar(p::PT, ::Type{T}) where {PT<:MP.AbstractPolynomial,T}
    return convert(MP.similar_type(PT, T), copy(p)) # Missing the `copy` in MP
end
Base.iterate(t::MP.Term) = iterate(t, 1)
function Base.iterate(t::MP.Term, state)
    if state == 1
        return MP.monomial(t), 2
    elseif state == 2
        return MP.coefficient(t), 3
    else
        return nothing
    end
end
function MA.operate!(
    ::SA.UnsafeAddMul{typeof(*)},
    mc::MP.AbstractPolynomial,
    val,
    c::MP.AbstractPolynomial,
)
    return MA.operate!(MA.add_mul, mc, val, c)
end
function MA.operate!(::typeof(SA.canonical), ::MP.AbstractPolynomial) end
# TODO Move to SA
function MA.promote_operation(
    ::typeof(SA.canonical),
    ::Type{P},
) where {P<:MP.AbstractPolynomialLike}
    return P
end

struct Polynomial{B,M<:MP.AbstractMonomial}
    monomial::M
    function Polynomial{B}(mono::MP.AbstractMonomial) where {B}
        return new{B,typeof(mono)}(mono)
    end
end

function Polynomial{B}(v::MP.AbstractVariable) where {B}
    return Polynomial{B}(MP.monomial(v))
end

function _algebra_element(p, ::Type{B}) where {B}
    basis = FullBasis{B,MP.monomial_type(p)}()
    return SA.AlgebraElement(
        p,
        SA.StarAlgebra(Polynomial{B}(MP.constant_monomial(p)), basis),
    )
end

function Base.:*(a::Polynomial{B}, b::Polynomial{B}) where {B}
    return _algebra_element(Mul{B}()(a.monomial, b.monomial), B)
end

function _show(io::IO, mime::MIME, p::Polynomial{B}) where {B}
    print(io, B)
    print(io, "(")
    show(io, mime, p.monomial)
    return print(io, ")")
end
function Base.show(io::IO, mime::MIME"text/plain", p::Polynomial)
    return _show(io, mime, p)
end
function Base.show(io::IO, p::Polynomial)
    return show(io, MIME"text/plain"(), p)
end

function Base.zero(::Type{Polynomial{B,M}}) where {B,M}
    return _algebra_element(zero(MP.polynomial_type(M, Rational{Int})), B)
end

Base.zero(p::Polynomial) = zero(typeof(p))

struct Mul{B} <: SA.MultiplicativeStructure end
