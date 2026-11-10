export LazyMatrix, macaulay, expand!

"""
    @enum(ShiftStatus, NOT_INCLUDED, INCLUDED, NOT_REDUNDANT, REDUNDANT)

* `NOT_INCLUDED` : Not included yet in the rows of the Macaulay matrix
* `INCLUDED` : Already included in the rows of the Macaulay matrix but redundancy unknown
* `NOT_REDUNDANT` : Already included in the rows of the Macaulay matrix and no redundancy with the other rows included
* `REDUNDANT` : Redundant with another row already included in the Macaulay matrix
"""
@enum(ShiftStatus, NOT_INCLUDED, INCLUDED, NOT_REDUNDANT, REDUNDANT)

"""
    mutable struct LazyMatrix{T,P,V,B}
        polynomials::V
        row_shifts::B
        shift_statuses::Vector{Vector{ShiftStatus}}
        column_basis::B
    end

Macaulay matrix of the polynomials `polynomials` in the variables
`MP.variables(column_basis)`. The rows are the
polynomials `shift * polynomials[j]` for each `shift` in `row_shifts`
and `j` such that `shift_statuses[i][j]` is `INCLUDED` or `NOT_REDUNDANT`
where `i` is the index of `shift` in `row_shifts`.
The columns are indexed by the monomials of `column_basis`.
The rows can be added with [`expand!`](@ref) and the matrix can be built with
`SparseArrays.sparse`.
"""
mutable struct LazyMatrix{
    T,
    P<:MP.AbstractPolynomialLike,
    V<:AbstractVector{P},
    B<:SubBasis{Monomial},
}
    polynomials::V
    row_shifts::B
    shift_statuses::Vector{Vector{ShiftStatus}}
    column_basis::B
end

function LazyMatrix(
    polynomials::V,
    variables = MP.variables(polynomials),
) where {T,P<:MP.AbstractPolynomialLike{T},V<:AbstractVector{P}}
    full = FullBasis{Monomial}(variables)
    empty = SA.SubBasis(full, SA.key_type(full)[])
    return LazyMatrix{T,P,V,typeof(empty)}(
        polynomials,
        empty,
        Vector{ShiftStatus}[],
        empty,
    )
end

function MP.monomial_type(::Type{<:LazyMatrix{T,P}}) where {T,P}
    return MP.monomial_type(P)
end

MP.variables(M::LazyMatrix) = MP.variables(M.column_basis)

function Base.show(io::IO, M::LazyMatrix)
    num_rows, num_cols = size(M)
    println(io, "$(num_rows)×$(num_cols) Macaulay matrix for polynomials:")
    for poly in M.polynomials
        println(io, "  ", poly)
    end
    println(io, "The row shifts are:")
    println(io, M.row_shifts)
    println(io, "The column basis is:")
    println(io, M.column_basis)
    return
end

function Base.size(M::LazyMatrix, i::Int)
    if i == 1
        return sum(M.shift_statuses; init = 0) do statuses
            return count(statuses) do s
                return s == INCLUDED || s == NOT_REDUNDANT
            end
        end
    elseif i == 2
        return length(M.column_basis)
    else
        return 1
    end
end

function Base.size(M::LazyMatrix)
    return (size(M, 1), size(M, 2))
end

function _sub_basis(M::LazyMatrix, monos)
    return parent(M.column_basis)[MP.monomial_vector(collect(monos))]
end

abstract type AbstractShiftsSelector end

"""
    struct ColumnDegrees{V} <: AbstractShiftsSelector
        degrees::V
    end

Select the shifts such that the maxdegree of the shifted polynomial is in
`degrees`.
"""
struct ColumnDegrees{V} <: AbstractShiftsSelector
    degrees::V
end

function select_shifts(vars, poly, d::ColumnDegrees)
    deg = MP.maxdegree(poly)
    if deg > maximum(d.degrees)
        return nothing
    else
        return MP.monomials(vars, max.(0, d.degrees .- deg))
    end
end

"""
    struct LeadingTargetColumns{M} <: AbstractShiftsSelector
        targets::Vector{M}
    end

Select the shifts such that the leading monomial of the shifted polynomial is
in `targets`.
"""
struct LeadingTargetColumns{M} <: AbstractShiftsSelector
    targets::Vector{M}
end

function select_shifts(vars, poly, d::LeadingTargetColumns{M}) where {M}
    mono = MP.leading_monomial(poly)
    return M[
        MP.div_multiple(target, mono) for
        target in d.targets if MP.divides(mono, target)
    ]
end

"""
    struct TargetColumns{M} <: AbstractShiftsSelector
        targets::Vector{M}
    end

Select the shifts such that one of the monomials of the shifted polynomial is
in `targets`.
"""
struct TargetColumns{M} <: AbstractShiftsSelector
    targets::Vector{M}
end

function select_shifts(vars, poly, d::TargetColumns{M}) where {M}
    return promote_type(M, MP.monomial_type(poly))[
        MP.div_multiple(target, mono) for mono in MP.monomials(poly) for
        target in d.targets if MP.divides(mono, target)
    ]
end

"""
    struct FixedShifts{M} <: AbstractShiftsSelector
        shifts::Vector{M}
    end

Select the shifts `shifts` for every polynomial.
"""
struct FixedShifts{M} <: AbstractShiftsSelector
    shifts::Vector{M}
end

function select_shifts(vars, poly, d::FixedShifts)
    return d.shifts
end

"""
    is_forever_trivial(M::LazyMatrix, col::MP.AbstractMonomialLike)

Return `true` if there does not exists any shift (even outside of
`M.row_shifts`) that such that the shift of any of the polynomial
has `col` as one of its monomial.
"""
function is_forever_trivial(M::LazyMatrix, col::MP.AbstractMonomialLike)
    return !any(M.polynomials) do poly
        return any(Base.Fix2(MP.divides, col), MP.monomials(poly))
    end
end

"""
    is_saturated(M::LazyMatrix, col::MP.AbstractMonomialLike)

Return `true` if all possible shift of a polynomial such that
`col` is one of the shifted monomial have been included.
"""
function is_saturated(M::LazyMatrix, col::MP.AbstractMonomialLike)
    for (j, poly) in enumerate(M.polynomials)
        for mono in MP.monomials(poly)
            if MP.divides(mono, col)
                shift = MP.div_multiple(col, mono)
                i = _key_index(M.row_shifts, shift)
                if isnothing(i) || M.shift_statuses[i][j] == NOT_INCLUDED
                    return false
                end
            end
        end
    end
    return true
end

"""
    expand!(M::LazyMatrix, shifts_selector::AbstractShiftsSelector; sparse_columns::Bool = true)

Add the row shifts selected by `shifts_selector`.
If `sparse_columns` is `false` then all monomials of degree up to the maximal
degree of the columns will be added to the columns. Otherwise, only the columns
that correspond to the monomial of one of the shifted polynomials will be added.
Return the number of rows added.
"""
function expand!(M::LazyMatrix, shifts_selector; sparse_columns::Bool = true)
    vars = MP.variables(M)
    MT = MP.monomial_type(typeof(M))
    row_monos_to_add = Dict{MT,Vector{ShiftStatus}}()
    if sparse_columns
        col_monos_to_add = Set{MT}()
    else
        col_maxdeg = isempty(M.column_basis) ? 0 : MP.maxdegree(M.column_basis)
    end
    added = 0
    for (j, poly) in enumerate(M.polynomials)
        selected_shifts = select_shifts(vars, poly, shifts_selector)
        if isnothing(selected_shifts)
            continue
        end
        for shift in selected_shifts
            i = _key_index(M.row_shifts, shift)
            if isnothing(i)
                if !haskey(row_monos_to_add, shift)
                    row_monos_to_add[shift] =
                        fill(NOT_INCLUDED, length(M.polynomials))
                end
                statuses = row_monos_to_add[shift]
            else
                statuses = M.shift_statuses[i]
            end
            if statuses[j] == NOT_INCLUDED
                added += 1
                statuses[j] = INCLUDED
                for mono in MP.monomials(poly)
                    col = shift * mono
                    if isnothing(_key_index(M.column_basis, col))
                        if sparse_columns
                            push!(col_monos_to_add, col)
                        else
                            col_maxdeg = max(col_maxdeg, MP.degree(col))
                        end
                    end
                end
            end
        end
    end
    if sparse_columns
        if !isempty(col_monos_to_add)
            M.column_basis =
                SA.merge_bases(M.column_basis, _sub_basis(M, col_monos_to_add))
        end
    else
        M.column_basis = maxdegree_basis(parent(M.column_basis), col_maxdeg)
    end
    if !isempty(row_monos_to_add)
        old_shifts = M.row_shifts
        M.row_shifts =
            SA.merge_bases(M.row_shifts, _sub_basis(M, keys(row_monos_to_add)))
        old_statuses = M.shift_statuses
        M.shift_statuses =
            Vector{Vector{ShiftStatus}}(undef, length(M.row_shifts))
        for (shift, statuses) in zip(old_shifts.keys, old_statuses)
            i = SA.key_index(M.row_shifts, shift)
            M.shift_statuses[i] = statuses
        end
        for (shift, statuses) in row_monos_to_add
            i = _key_index(M.row_shifts, shift)
            M.shift_statuses[i] = statuses
        end
    end
    return added
end

"""
    macaulay(polynomials, maxdegree, variables = MP.variables(polynomials); kws...)

Return the [`LazyMatrix`](@ref) of the Macaulay matrix of degree `maxdegree`
of `polynomials` in the variables `variables`.
The keyword arguments are passed to [`expand!`](@ref).
"""
function macaulay(
    polynomials,
    maxdegree,
    variables = MP.variables(polynomials);
    kws...,
)
    M = LazyMatrix(polynomials, variables)
    expand!(M, ColumnDegrees(0:maxdegree); kws...)
    return M
end

function SparseArrays.sparse(M::LazyMatrix{T}) where {T}
    row = 0
    I = Int[]
    J = Int[]
    K = T[]
    for (shift, statuses) in
        zip(keys_as_monomials(M.row_shifts), M.shift_statuses)
        for (j, status) in enumerate(statuses)
            if status == INCLUDED || status == NOT_REDUNDANT
                row += 1
                for t in MP.terms(M.polynomials[j])
                    push!(I, row)
                    push!(J, _key_index(M.column_basis, shift * MP.monomial(t)))
                    push!(K, MP.coefficient(t))
                end
            end
        end
    end
    @assert row == size(M, 1)
    return SparseArrays.sparse(I, J, K, row, size(M, 2))
end
