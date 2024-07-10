"""
    abstract type AbstractMultipleOrthogonal <: AbstractMonomialIndexed end

Polynomial basis such that ``\\langle p_i(x), p_j(x) \\rangle = 0`` if ``i \\neq j`` where
```math
\\langle p(x), q(x) \\rangle = \\int p(x)q(x) w(x) dx
```
where the weight is a product of weight functions
``w(x) = w_1(x_1)w_2(x_2) \\cdots w_n(x_n)`` in each variable.
The polynomial of the basis are product of univariate polynomials:
``p(x) = p_1(x_1)p_2(x_2) \\cdots p_n(x_n)``.
where the univariate polynomials of variable `x_i` form an univariate
orthogonal basis for the weight function `w_i(x_i)`.
Therefore, they satisfy the recurrence relation
```math
d_k p_k(x_i) = (a_k x_i + b_k) p_{k-1}(x_i) + c_k p_{k-2}(x_i)
```
where [`reccurence_first_coef`](@ref) gives `a_k`,
[`reccurence_second_coef`](@ref) gives `b_k`,
[`reccurence_third_coef`](@ref) gives `c_k` and
[`reccurence_deno_coef`](@ref) gives `d_k`.
"""
abstract type AbstractMultipleOrthogonal <: AbstractMonomialIndexed end

"""
    reccurence_first_coef(B::Type{<:AbstractMultipleOrthogonal}, degree::Integer)

Return `a_{degree}` in recurrence equation
```math
d_k p_k(x_i) = (a_k x_i + b_k) p_{k-1}(x_i) + c_k p_{k-2}(x_i)
```
"""
function reccurence_first_coef end

"""
    reccurence_second_coef(B::Type{<:AbstractMultipleOrthogonal}, degree::Integer)

Return `b_{degree}` in recurrence equation
```math
d_k p_k(x_i) = (a_k x_i + b_k) p_{k-1}(x_i) + c_k p_{k-2}(x_i)
```
"""
function reccurence_second_coef end

"""
    reccurence_third_coef(B::Type{<:AbstractMultipleOrthogonal}, degree::Integer)

Return `c_{degree}` in recurrence equation
```math
d_k p_k(x_i) = (a_k x_i + b_k) p_{k-1}(x_i) + c_k p_{k-2}(x_i)
```
"""
function reccurence_third_coef end

"""
    reccurence_deno_coef(B::Type{<:AbstractMultipleOrthogonal}, degree::Integer)

Return `d_{degree}` in recurrence equation
```math
d_k p_k(x_i) = (a_k x_i + b_k) p_{k-1}(x_i) + c_k p_{k-2}(x_i)
```
"""
function reccurence_deno_coef end

function univariate_eval!(
    values::Vector,
    value,
    ::Type{B},
) where {B}
    if 1 in eachindex(values)
        values[1] = one(value)
    end
    if 2 in eachindex(values)
        values[1] = one(value)
    end
    for d in 2:(length(values) - 1)
    end
end

"""
    univariate_orthogonal_basis(
        B::Type{<:AbstractMultipleOrthogonal},
        variable::MP.AbstractVariable,
        degree::Integer,
    )

Return the vector of univariate polynomials of the basis `B` up to `degree`
with variable `variable`.
"""
function univariate_orthogonal_basis(
    B::Type{<:AbstractMultipleOrthogonal},
    variable::MP.AbstractVariable,
    degree::Integer,
)
    @assert degree >= 0
    if degree == 0
        return MP.polynomial_type(
            Polynomial{B,MP.monomial_type(variable)},
            Int,
        )[one(variable)]
    elseif degree == 1
        return push!(
            univariate_orthogonal_basis(B, variable, 0),
            degree_one_univariate_polynomial(B, variable),
        )
    else
        previous = univariate_orthogonal_basis(B, variable, degree - 1)
        a = reccurence_first_coef(B, degree)
        b = reccurence_second_coef(B, degree)
        c = reccurence_third_coef(B, degree)
        d = reccurence_deno_coef(B, degree)
        next = MA.@rewrite(
            (a * variable + b) * previous[degree] + c * previous[degree-1]
        )
        if !isone(d)
            next = next / d
        end
        push!(previous, next)
        return previous
    end
end

function _covering(::FullBasis{B,M}, monos) where {B,M}
    to_add = collect(monos)
    m = Set{M}(to_add)
    while !isempty(to_add)
        mono = pop!(to_add)
        for v in MP.variables(mono)
            step = even_odd_separated(B) ? 2 : 1
            vstep = v^step
            if MP.divides(vstep, mono)
                new_mono = MP.map_exponents(-, mono, vstep)
                if !(new_mono in m)
                    push!(m, new_mono)
                    push!(to_add, new_mono)
                end
            end
        end
    end
    return collect(m)
end

function explicit_basis_covering(
    full::FullBasis{BM,M},
    monos::SubBasis{B,M},
) where {BM<:AbstractMonomial,B<:AbstractMultipleOrthogonal,M}
    full = FullBasis{B,M}()
    return SubBasis{BM}(_covering(full, monos.monomials))
end

function explicit_basis_covering(
    full::FullBasis{B,M},
    monos::SubBasis{<:AbstractMonomial,M},
) where {B<:AbstractMultipleOrthogonal,M}
    return SubBasis{B}(_covering(full, monos.monomials))
end

function _scalar_product_function(::Type{<:AbstractMultipleOrthogonal}, i::Int) end

function LinearAlgebra.dot(p, q, basis_type::Type{<:AbstractMultipleOrthogonal})
    return _integral(p * q, basis_type)
end

function _integral(p::Number, basis_type::Type{<:AbstractMultipleOrthogonal})
    return p * _scalar_product_function(basis_type, 0)
end

function _integral(
    ::MP.AbstractVariable,
    basis_type::Type{<:AbstractMultipleOrthogonal},
)
    return _scalar_product_function(basis_type, 1)
end

function _integral(
    p::MP.AbstractMonomial,
    basis_type::Type{<:AbstractMultipleOrthogonal},
)
    return prod([
        _scalar_product_function(basis_type, i) for i in MP.exponents(p)
    ])
end

function _integral(
    p::MP.AbstractTerm,
    basis_type::Type{<:AbstractMultipleOrthogonal},
)
    return MP.coefficient(p) * _integral(MP.monomial(p), basis_type)
end

function _integral(
    p::MP.AbstractPolynomial,
    basis_type::Type{<:AbstractMultipleOrthogonal},
)
    return sum([_integral(t, basis_type) for t in MP.terms(p)])
end

function MP.coefficients(
    p,
    basis::SubBasis{B,M},
) where {B<:AbstractMultipleOrthogonal,M}
    poly_p = MP.polynomial(p)
    return map(basis) do el
        q = SA.coeffs(el, FullBasis{Monomial,M}())
        poly_q = MP.polynomial(q)
        return LinearAlgebra.dot(poly_p, poly_q, B) /
               LinearAlgebra.dot(poly_q, poly_q, B)
    end
end

function SA.coeffs(
    p::Polynomial{B},
    basis::SubBasis{Monomial},
) where {B<:AbstractMultipleOrthogonal}
    full = implicit_basis(basis)
    return SA.coeffs(algebra_element(SA.coeffs(p, full), full), basis)
end

function SA.coeffs(
    p::Polynomial{B,M},
    ::FullBasis{Monomial},
) where {B<:AbstractMultipleOrthogonal,M}
    return sparse_coefficients(
        prod(
            MP.powers(p.monomial);
            init = MP.constant_monomial(M),
        ) do (var, deg)
            return univariate_orthogonal_basis(B, var, deg)[deg+1]
        end,
    )
end
