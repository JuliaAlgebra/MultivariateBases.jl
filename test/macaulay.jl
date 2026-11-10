using Test
using DynamicPolynomials
import MultivariateBases as MB
import SparseArrays

function _test_rows(M::MB.LazyMatrix, expected_rows)
    A = Matrix(SparseArrays.sparse(M))
    @test size(A) == size(M)
    columns = MB.keys_as_monomials(M.column_basis)
    rows = [polynomial(A[i, :], columns) for i in axes(A, 1)]
    @test length(rows) == length(expected_rows)
    @test Set(rows) == Set(expected_rows)
end

@testset "Macaulay matrix" begin
    @polyvar x y
    p = x * y - 1
    q = x^2 - y
    M = MB.macaulay([p, q], 2)
    @test size(M) == (2, 4)
    @test size(M, 3) == 1
    _test_rows(M, [p, q])
    @test MB.keys_as_monomials(M.column_basis) == [1, y, x * y, x^2]
    @test MB.is_saturated(M, constant_monomial(p))
    @test !MB.is_saturated(M, x * y)
    @test !MB.is_forever_trivial(M, x^2 * y)
    @test MB.expand!(M, MB.ColumnDegrees(0:3)) == 4
    @test size(M) == (6, 9)
    _test_rows(M, [p, x * p, y * p, q, x * q, y * q])
    @test MB.is_saturated(M, y)
    @test MB.expand!(M, MB.ColumnDegrees(0:3)) == 0
    @test sprint(show, M) isa String
    M = MB.macaulay([p, q], 3; sparse_columns = false)
    @test size(M) == (6, 10)
    @test MB.keys_as_monomials(M.column_basis) == monomials([x, y], 0:3)
    @test isnothing(MB.select_shifts([x, y], p, MB.ColumnDegrees(0:1)))
    M = MB.macaulay([p, x^3 - 1], 2)
    _test_rows(M, [p])
    @polyvar z
    M = MB.macaulay([q], 3, [x, y, z])
    _test_rows(M, [q, x * q, y * q, z * q])
end

@testset "Shift selectors" begin
    @polyvar x y z
    p = x * y - 1
    q = x^2 - y
    @test !MB.is_forever_trivial(MB.LazyMatrix([p, q]), z)
    @test MB.is_forever_trivial(MB.LazyMatrix([q]), z)
    M = MB.LazyMatrix([p, q])
    @test MB.expand!(M, MB.LeadingTargetColumns([x^2 * y])) == 2
    _test_rows(M, [x * p, y * q])
    M = MB.LazyMatrix([p, q])
    @test MB.expand!(M, MB.TargetColumns([x^2 * y])) == 4
    _test_rows(M, [x * p, x^2 * y * p, y * q, x^2 * q])
    M = MB.LazyMatrix([p, q])
    @test MB.expand!(M, MB.FixedShifts([x])) == 2
    _test_rows(M, [x * p, x * q])
end
