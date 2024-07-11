"""
    struct Trigonometric <: AbstractMonomialIndexed end

Univariate trigonometric basis is
```
a0 + a1 cos(ωt) + a2 sin(ωt) + a3 cos(2ωt) + a4 sin(2ωt)
```
"""
struct Trigonometric <: AbstractMonomialIndexed end

_cos_id(d) = iszero(d) ? 0 : 2d - 1
_sin_id(d) = 2d
_is_cos(d) = isodd(d)
_is_sin(d) = isodd(d)
_id(d) = div(d + 1, 2)

# https://en.wikipedia.org/wiki/Chebyshev_polynomials#Properties
# Using
# sin(a + b) = sin(a) cos(b) + cos(a) sin(b)
# sin(a - b) = sin(a) cos(b) - cos(a) sin(b)
# If a > b
# sin(a) cos(b) = sin(a + b) + sin(a - b)
# sin(a) cos(b) = sin(a + b) + sin(a - b)
function univariate_mul!(::Mul{Trigonometric}, terms, var, a, b)
    @assert !iszero(a)
    @assert !iszero(b)
    I = eachindex(terms)
    da = _id(a)
    db = _id(b)
    for i in I
        if _is_cos(a) == _is_cos(b)
            # Chebyshev first kind
            mono = MP.monomial(terms[i]) * var^(_cos_id(da + db))
            terms[i] = MA.mul!!(terms[i], var^_cos_id(abs(da - db)))
            terms[i] = MA.operate!!(/, terms[i], 2)
            α = MA.copy_if_mutable(MP.coefficient(terms[i]))
            push!(terms, MP.term(α, mono))
            # cos(a + b) = cos(a) cos(b) - sin(a) sin(b)
            # cos(a - b) = cos(a) cos(b) + sin(a) sin(b)
            # cos(a - b) - cos(a + b)
            if _is_sin(a)
                terms[end] = MA.operate!!(*, terms[end], -1)
            end
        else
            if _is_cos(a)
                da, db = db, da
            end
            # sin(da) * cos(db)
            if da == db
                # sin(da) * cos(da) = sin(2da) / 2
                terms[i] = MA.mul!!(terms[i], var^_cos_id(da + db))
                terms[i] = MA.operate!!(/, terms[i], 2)
            else
                # Using
                # sin(a + b) = sin(a) cos(b) + cos(a) sin(b)
                # sin(a - b) = sin(a) cos(b) - cos(a) sin(b)
                # If a > b
                # sin(a) cos(b) = (sin(a + b) + sin(a - b)) / 2
                # If a < b
                # sin(a) cos(b) = (sin(b + a) - sin(b - a)) / 2
                mono = MP.monomial(terms[i]) * var^(_sin_id(da + db))
                terms[i] = MA.mul!!(terms[i], var^_sin_id(abs(da - db)))
                terms[i] = MA.operate!!(/, terms[i], 2)
                α = MA.copy_if_mutable(MP.coefficient(terms[i]))
                push!(terms, MP.term(α, mono))
                if da < db
                    terms[i] = MA.operate!!(*, terms[i], -1)
                end
            end
        end
    end
    return
end

function degree_one_univariate_polynomial(::Type{Trigonometric}, variable)
    MA.@rewrite(variable + 0)
end

function recurrence_eval(
    ::Type{Trigonometric},
    previous::AbstractVector,
    value,
    degree,
)
    d = _id(degree)
    if _is_cos(degree)
        # Chebyshev first order
        return 2 * value * previous[_cos_id(d - 1) + 1] - previous[_cos_id(d - 2) + 1]
    else
        return sqrt(1 - previous[degree]^2)
    end
end

function _promote_coef(::Type{T}, ::Type{Trigonometric}) where {T}
    return _promote_div(T)
end
