# Basis-specific operations that currently expand to the monomial basis.
const _BasisPolynomial{T} = SA.AlgebraElement{
    T,
    <:SA.StarAlgebra{
        <:Variables{<:Union{ScaledMonomial,AbstractMultipleOrthogonal}},
    },
}

function _monomial_expansion(p::_BasisPolynomial)
    return convert_basis(FullBasis{Monomial}(MP.variables(p)), p)
end

const _ScaledPolynomial{T} =
    SA.AlgebraElement{T,<:SA.StarAlgebra{<:Variables{ScaledMonomial}}}

# The conversion between the scaled monomial and the monomial basis is a
# diagonal scaling so the substitution is computed in the monomial basis and
# the result is converted back to the basis family of `p`.
_same_family(_, value) = value

function _same_family(
    ::SA.StarAlgebra{<:Variables{B}},
    q::MP.AbstractPolynomialLike,
) where {B}
    return convert_basis(FullBasis{B}(MP.variables(q)), q)
end

function MP.substitute(
    s::MP.AbstractSubstitutionType,
    p::_ScaledPolynomial,
    args::MP.Substitutions,
)
    return _same_family(
        parent(p),
        MP.substitute(s, _monomial_expansion(p), args),
    )
end

# The conversion between an orthogonal basis and the monomial basis is
# ill-conditioned so the substitution is computed in the orthogonal basis using
# the three-term recurrence, the substituted values are only multiplied
# through `_mul_variable`.
const _OrthogonalPolynomial{T} = SA.AlgebraElement{
    T,
    <:SA.StarAlgebra{<:Variables{<:AbstractMultipleOrthogonal}},
}

function _push_scaled!(res, q, num, den)
    for (k, v) in SA.nonzero_pairs(q)
        SA.unsafe_push!(res, k, v * num / den)
    end
    return res
end

# Coefficients of `x_i * q` where `q` contains coefficients in the basis `B`
function _mul_variable(::Type{B}, q, i) where {B}
    res = zero(q)
    # `p_1(x) = α x + β`
    β = degree_one_univariate_polynomial(B, 0)
    α = degree_one_univariate_polynomial(B, 1) - β
    for (exp, v) in SA.nonzero_pairs(q)
        n = exp[i]
        if iszero(n)
            # `x p_0 = (p_1 - β p_0) / α`
            SA.unsafe_push!(res, _increment(exp, 1, i), v / α)
            if !iszero(β)
                SA.unsafe_push!(res, exp, -v * β / α)
            end
        else
            # `x p_n = (d p_{n+1} - b p_n - c p_{n-1}) / a`
            a = reccurence_first_coef(B, n + 1)
            b = reccurence_second_coef(B, n + 1)
            c = reccurence_third_coef(B, n + 1)
            d = reccurence_deno_coef(B, n + 1)
            SA.unsafe_push!(res, _increment(exp, 1, i), v * d / a)
            if !iszero(b)
                SA.unsafe_push!(res, exp, -v * b / a)
            end
            SA.unsafe_push!(res, _increment(exp, -1, i), -v * c / a)
        end
    end
    return MA.operate!(SA.canonical, res)
end

# Coefficients of `value * q` where `value` is in the monomial basis
_mul_value(::Type, q, value::Number, _) = _push_scaled!(zero(q), q, value, 1)

function _mul_value(::Type{B}, q, value, vars) where {B}
    res = zero(q)
    for t in MP.terms(value)
        r = q
        for (var, deg) in MP.powers(MP.monomial(t))
            i = findfirst(isequal(var), vars)
            for _ in 1:deg
                r = _mul_variable(B, r, i)
            end
        end
        _push_scaled!(res, r, MP.coefficient(t), 1)
    end
    return MA.operate!(SA.canonical, res)
end

# Coefficients of `p_n(value) * q` using the three-term recurrence
# `d p_k = (a value + b) p_{k-1} + c p_{k-2}`
function _mul_univariate(::Type{B}, q, n, value, vars) where {B}
    if iszero(n)
        return q
    end
    β = degree_one_univariate_polynomial(B, 0)
    α = degree_one_univariate_polynomial(B, 1) - β
    prev = q
    cur = _push_scaled!(zero(q), _mul_value(B, q, value, vars), α, 1)
    cur = MA.operate!(SA.canonical, _push_scaled!(cur, q, β, 1))
    for k in 2:n
        a = reccurence_first_coef(B, k)
        b = reccurence_second_coef(B, k)
        c = reccurence_third_coef(B, k)
        d = reccurence_deno_coef(B, k)
        next = _push_scaled!(zero(q), _mul_value(B, cur, value, vars), a, d)
        _push_scaled!(next, cur, b, d)
        _push_scaled!(next, prev, c, d)
        prev, cur = cur, MA.operate!(SA.canonical, next)
    end
    return cur
end

_monomial_value(value::Number) = value

_coefficient_type(value::Number) = typeof(value)
_coefficient_type(value) = MP.coefficient_type(value)

function _monomial_value(value)
    q = MP.polynomial(value)
    return convert_basis(FullBasis{Monomial}(MP.variables(q)), q)
end

function MP.substitute(
    st::MP.AbstractSubstitutionType,
    p::_OrthogonalPolynomial,
    args::MP.Substitutions,
)
    foreach(MP._check_substitution, args)
    B = _family(parent(p))
    values = [
        _monomial_value(MP._substitution_value(Any, st, v, args...)) for
        v in MP.variables(p)
    ]
    polys = filter(v -> !(v isa Number), values)
    one_mono = prod(
        MP.constant_monomial,
        polys;
        init = MP.constant_monomial(MP.monomial_type(p)),
    )
    vars = MP.variables(one_mono)
    full = FullBasis{B}(vars)
    T = _promote_coef(
        MA.promote_operation(
            *,
            MP.coefficient_type(p),
            promote_type(map(_coefficient_type, values)...),
        ),
        B,
    )
    one_key = MP.exponents(one_mono)
    res = SA.SparseCoefficients(typeof(one_key)[], T[], SA.comparable(full))
    for (exp, v) in SA.nonzero_pairs(SA.coeffs(p))
        q = SA.SparseCoefficients([one_key], T[v], SA.comparable(full))
        for (n, value) in zip(exp, values)
            q = _mul_univariate(B, q, n, value, vars)
        end
        _push_scaled!(res, q, 1, 1)
    end
    MA.operate!(SA.canonical, res)
    if st isa MP.Eval && isempty(polys)
        return isempty(SA.keys(res)) ? zero(T) : only(SA.values(res))
    end
    return algebra_element(res, full)
end

_family(::SA.StarAlgebra{<:Variables{B}}) where {B} = B

function MP.subs(p::_BasisPolynomial, args::MP.AbstractSubstitution...)
    return MP.substitute(MP.Subs(), p, args)
end

function (p::_BasisPolynomial)(args::MP.AbstractSubstitution...)
    return MP.substitute(MP.Eval(), p, args)
end

function (p::_BasisPolynomial)(x::NTuple{N,<:Number}) where {N}
    return p(collect(x))
end

function (p::_BasisPolynomial)(x::AbstractVector{<:Number})
    return p(MP.variables(p) => x)
end

function (p::_BasisPolynomial)(x::Number...)
    return p(x)
end

function MP.differentiate(p::_BasisPolynomial, x::MP.AbstractVariable)
    return MP.differentiate(_monomial_expansion(p), x)
end
