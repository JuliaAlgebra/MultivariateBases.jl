"""
    abstract type AbstractMultipleOrthogonalBasis{P} <: AbstractPolynomialVectorBasis{P, Vector{P}} end

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
abstract type AbstractMultipleOrthogonalBasis{P} <: AbstractPolynomialVectorBasis{P, Vector{P}} end

"""
    reccurence_first_coef(B::Type{<:AbstractMultipleOrthogonalBasis}, degree::Integer)

Return `a_{degree}` in recurrence equation
```math
d_k p_k(x_i) = (a_k x_i + b_k) p_{k-1}(x_i) + c_k p_{k-2}(x_i)
```
"""
function reccurence_first_coef end

"""
    reccurence_second_coef(B::Type{<:AbstractMultipleOrthogonalBasis}, degree::Integer)

Return `b_{degree}` in recurrence equation
```math
d_k p_k(x_i) = (a_k x_i + b_k) p_{k-1}(x_i) + c_k p_{k-2}(x_i)
```
"""
function reccurence_second_coef end

"""
    reccurence_third_coef(B::Type{<:AbstractMultipleOrthogonalBasis}, degree::Integer)

Return `c_{degree}` in recurrence equation
```math
d_k p_k(x_i) = (a_k x_i + b_k) p_{k-1}(x_i) + c_k p_{k-2}(x_i)
```
"""
function reccurence_third_coef end

"""
    reccurence_deno_coef(B::Type{<:AbstractMultipleOrthogonalBasis}, degree::Integer)

Return `d_{degree}` in recurrence equation
```math
d_k p_k(x_i) = (a_k x_i + b_k) p_{k-1}(x_i) + c_k p_{k-2}(x_i)
```
"""
function reccurence_deno_coef end

"""
    univariate_orthogonal_basis(B::Type{<:AbstractMultipleOrthogonalBasis},
                                variable::MP.AbstractVariable, degree::Integer)

Return the vector of univariate polynomials of the basis `B` up to `degree`
with variable `variable`.
"""
function univariate_orthogonal_basis(
    B::Type{<:AbstractMultipleOrthogonalBasis}, variable::MP.AbstractVariable, degree::Integer)

    @assert degree >= 0
    if degree == 0
        return polynomial_type(B, typeof(variable))[one(variable)]
    elseif degree == 1
        return push!(univariate_orthogonal_basis(B, variable, 0),
                     degree_one_univariate_polynomial(B, variable))
    else
        previous = univariate_orthogonal_basis(B, variable, degree - 1)
        a = reccurence_first_coef(B, degree)
        b = reccurence_second_coef(B, degree)
        c = reccurence_third_coef(B, degree)
        d = reccurence_deno_coef(B, degree)
        next = MA.@rewrite((a * variable + b) * previous[degree] + c * previous[degree - 1])
        if !isone(d)
            next = next / d
        end
        push!(previous, next)
        return previous
    end
end

function _basis_from_monomials(B::Type{<:AbstractMultipleOrthogonalBasis}, variables, monos)
    univariate = [univariate_orthogonal_basis(
                      B, variable, maximum(mono -> degree(mono, variable), monos))
                  for variable in variables]
    return B([
        prod(i -> univariate[i][degree(mono, variables[i]) + 1],
             eachindex(variables))
        for mono in monos])
end

function maxdegree_basis(B::Type{<:AbstractMultipleOrthogonalBasis}, variables, maxdegree::Int)
    return _basis_from_monomials(B, variables, MP.monomials(variables, 0:maxdegree))
end

function basis_covering_monomials(B::Type{<:AbstractMultipleOrthogonalBasis}, monos::AbstractVector{<:AbstractMonomial})
    to_add = collect(monos)
    m = Set(monos)
    while !isempty(to_add)
        mono = pop!(to_add)
        for v in MP.variables(mono)
            step = even_odd_separated(B) ? 2 : 1
            vstep = v^step
            if MP.divides(vstep, mono)
                new_mono = MP.mapexponents(-, mono, vstep)
                if !(new_mono in m)
                    push!(m, new_mono)
                    push!(to_add, new_mono)
                end
            end
        end
    end
    return _basis_from_monomials(B, variables(monos), MP.monovec(collect(m)))
end

function scalar_product_function(::Type{<:AbstractMultipleOrthogonalBasis}, i::Int) end

LinearAlgebra.dot(p, q, basis_type::Type{<:AbstractMultipleOrthogonalBasis}) = _integral(p*q, basis_type)

function _integral(p::Number, basis_type::Type{<:AbstractMultipleOrthogonalBasis})
    return p*scalar_product_function(basis_type, 0)
end

function _integral(p::MP.AbstractVariable, basis_type::Type{<:AbstractMultipleOrthogonalBasis})
    return scalar_product_function(basis_type, 1)
end

function _integral(p::MP.AbstractMonomial, basis_type::Type{<:AbstractMultipleOrthogonalBasis})
    return prod([scalar_product_function(basis_type, i) for i in exponents(p)])
end

function _integral(p::MP.AbstractTerm, basis_type::Type{<:AbstractMultipleOrthogonalBasis})
    return coefficient(p)*_integral(monomial(p), basis_type)
end

function _integral(p::MP.AbstractPolynomial, basis_type::Type{<:AbstractMultipleOrthogonalBasis})
    return sum([_integral(t, basis_type) for t in terms(p)])
end

function MP.coefficients(p, basis::AbstractMultipleOrthogonalBasis; check = true)
    B = typeof(basis)
    coeffs = [LinearAlgebra.dot(p, el, B)/LinearAlgebra.dot(el, el, B) for el in basis]
    idx = findall(c -> !isapprox(c, 0, atol = 1e-10), coeffs)
    return coeffs[idx]
end 
