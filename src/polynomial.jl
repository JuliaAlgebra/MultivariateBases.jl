# TODO Add to MultivariatePolynomials
MP.variables(p::SA.AlgebraElement) = MP.variables(explicit_basis(p))
Base.keytype(p::MP.AbstractPolynomialLike) = MP.monomial_type(p)
SA.value_type(p::MP.AbstractPolynomialLike) = MP.coefficient_type(p)
#Base.keys(p::MP.AbstractPolynomial) = MP.monomials(p)
SA.nonzero_pairs(p::MP.AbstractPolynomialLike) = MP.terms(p)
function Base.similar(p::PT, ::Type{T}) where {PT<:MP.AbstractPolynomial,T}
    return convert(MP.similar_type(PT, T), copy(p)) # Missing the `copy` in MP
end
function Base.getindex(p::MP.AbstractPolynomialLike, mono::MP.AbstractMonomial)
    return MP.coefficient(p, mono)
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
function SA.unsafe_push!(p::MP.AbstractPolynomial, mono::MP.AbstractMonomial, α)
    return MA.operate!(MA.add_mul, p, α, mono)
end
function MA.operate!(
    ::SA.UnsafeAddMul{typeof(*)},
    mc::MP.AbstractPolynomial,
    val,
    c::MP.AbstractPolynomialLike,
)
    return MA.operate!(MA.add_mul, mc, val, c)
end
MA.operate!(::typeof(SA.canonical), p::MP.AbstractPolynomial) = p
function MA.promote_operation(
    ::typeof(SA.canonical),
    ::Type{P},
) where {P<:MP.AbstractPolynomialLike}
    return P
end

abstract type AbstractMonomialIndexed end

"""
    struct Polynomial{B<:AbstractMonomialIndexed,M<:MP.AbstractMonomial}
        monomial::M
        function Polynomial{B}(mono::MP.AbstractMonomial) where {B}
            return new{B,typeof(mono)}(mono)
        end
    end

Polynomial of basis `FullBasis{B,M}()` at index `monomial`.
"""
struct Polynomial{B<:AbstractMonomialIndexed,V,E}
    variables::Variables{B,V}
    exponents::E
end

function Polynomial(v::Variables{B,V}, e) where {B,V}
    return Polynomial{B,V,typeof(e)}(v, e)
end

function Polynomial{B}(v::MP.AbstractVariable) where {B}
    vars = MP.variables(v)
    # MP.exponents(v) gives (1,) for DynamicPolynomials so we use
    # a custom `map` as workaround
    return Polynomial(Variables{B}(vars), map(_ -> 1, vars))
end

MP.exponents(p::Polynomial) = p.exponents

MP.monomial(p::Polynomial) = MP.monomial(MP.variables(p), MP.exponents(p))

function Base.hash(p::Polynomial{B}, u::UInt) where {B}
    return hash(p.variables, hash(p.monomial, u))
end

function Base.isequal(p::Polynomial{B}, q::Polynomial{B}) where {B}
    return isequal(p.variables, q.variables) && isequal(p.exponents, q.exponents)
end

Base.isone(p::Polynomial) = all(iszero, p.exponents)

# Needed for `BoundsError`
Base.iterate(p::Polynomial) = p, nothing
Base.iterate(::Polynomial, ::Nothing) = nothing

function Base.:(==)(p::Polynomial{B}, q::Polynomial{B}) where {B}
    return p.variables == q.variables && p.exponents == q.exponents
end

MP.variables(p::Polynomial) = MP.variables(p.variables)
MP.nvariables(p::Polynomial) = MP.nvariables(p.variables)

MP.monomial_type(::Type{<:SA.SparseCoefficients{K}}) where {K} = K
MP.polynomial(p::Polynomial) = MP.polynomial(algebra_element(p))

function algebra_element(p, basis::SA.AbstractBasis)
    return SA.AlgebraElement(p, algebra(basis))
end

function _algebra_element(p, ::Type{B}) where {B<:AbstractMonomialIndexed}
    return algebra_element(
        sparse_coefficients(p),
        FullBasis{B}(MP.variables(p)),
    )
end

function algebra_element(p::Polynomial{B}) where {B}
    return _algebra_element(MP.monomial(p), B)
end

function Base.:*(a::Polynomial{B}, b::SA.AlgebraElement) where {B}
    return _algebra_element(a) * b
end

function _show(io::IO, mime::MIME, p::Polynomial{B}) where {B}
    if B != Monomial
        print(io, B)
        print(io, "(")
    end
    print(io, SA.trim_LaTeX(mime, sprint(show, mime, MP.monomial(p.variables.variables, p.exponents))))
    if B != Monomial
        print(io, ")")
    end
    return
end

function Base.show(io::IO, mime::MIME"text/latex", p::Polynomial)
    print(io, "\$\$ ")
    _show(io, mime, p)
    print(io, " \$\$")
    return
end

function Base.show(io::IO, mime::MIME"text/plain", p::Polynomial)
    return _show(io, mime, p)
end

function Base.show(io::IO, mime::MIME"text/print", p::Polynomial)
    return _show(io, mime, p)
end

Base.show(io::IO, p::Polynomial) = show(io, MIME"text/plain"(), p)
Base.print(io::IO, p::Polynomial) = show(io, MIME"text/print"(), p)

function Base.zero(::Type{Polynomial{B,M}}) where {B,M}
    return _algebra_element(zero(MP.polynomial_type(M, Rational{Int})), B)
end

Base.zero(p::Polynomial) = zero(typeof(p))

function convert_basis(basis::SA.AbstractBasis, p::MP.AbstractPolynomialLike)
    return convert_basis(basis, _algebra_element(p, Monomial))
end

function convert_basis(basis::SA.AbstractBasis, p::SA.AlgebraElement)
    return SA.AlgebraElement(SA.coeffs(p, basis), algebra(basis))
end

function Base.isapprox(
    p::MP.AbstractPolynomialLike,
    a::SA.AlgebraElement;
    kws...,
)
    return isapprox(p, MP.polynomial(a); kws...)
end

function Base.isapprox(a::SA.AlgebraElement, b::SA.AlgebraElement; kws...)
    return isapprox(MP.polynomial(a), b; kws...)
end

function Base.isapprox(
    a::SA.AlgebraElement,
    p::MP.AbstractPolynomialLike;
    kws...,
)
    return isapprox(p, a; kws...)
end

function Base.isapprox(a::SA.AlgebraElement, α::Number; kws...)
    return isapprox(
        a,
        α * constant_algebra_element(typeof(SA.basis(a)), typeof(α));
        kws...,
    )
end
