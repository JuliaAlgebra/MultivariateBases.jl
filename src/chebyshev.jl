struct ChebyshevBasis{P} <: AbstractPolynomialVectorBasis{P, Vector{P}}
    polynomials::Vector{P}
end

function chebyshev_polynomial_first_kind(variable::MP.AbstractVariable, degree::Integer)
    @assert degree >= 0
    if degree == 0
        return [MA.@rewrite(0 * variable + 1)]
    elseif degree == 1
        return push!(chebyshev_polynomial_first_kind(variable, 0),
                     MA.@rewrite(1 * variable + 0))
    else
        previous = chebyshev_polynomial_first_kind(variable, degree - 1)
        next = MA.@rewrite(2variable * previous[degree] - previous[degree - 1])
        push!(previous, next)
        return previous
    end
end

function maxdegree_basis(B::Type{ChebyshevBasis}, variables, maxdegree::Int)
    univariate = [chebyshev_polynomial_first_kind(variable, maxdegree) for variable in variables]
    return ChebyshevBasis([
        prod(i -> univariate[i][degree(mono, variables[i]) + 1],
             eachindex(variables))
        for mono in MP.monomials(variables, 0:maxdegree)])
end
