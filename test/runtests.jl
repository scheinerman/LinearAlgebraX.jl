using Test
using LinearAlgebra, LinearAlgebraX, Mods


@testset "Basics" begin
    A = big.(rand(Int, 10, 10) .% 100)

    include("../extras/hilbert.jl")

    H = hilbert(12)
    @test rankx(H) == 12

    @test detx(A) == det(Rational.(A))
    @test detx(A) == detx(A')
    @test A * invx(A) == eye(BigInt, 10)

    A = ones(Int, 3, 5)
    N = nullspacex(A)
    @test all(0 .== A * N)
end

@testset "Polynomials" begin
    using SimplePolynomials

    A = triu(ones(Int, 5, 5))
    p = char_poly(A)
    x = getx()
    @test p == (x - 1)^5

    A = ones(Int, 5, 5)
    p = char_poly(A)
    @test p == x^4 * (x - 5)
end

@testset "Modular" begin
    T = Mod{17}
    A = rand(T, 6, 6)
    p = char_poly(A)
    @test detx(A) == p(0)


    A = rand(Int, 20, 20) .% 100
    @test rankx(A) <= 20
end

@testset "Complex" begin
    A = rand(Int, 5, 5) .% 100
    B = rand(Int, 5, 5) .% 100
    M = A + im * B
    @test detx(M) == cofactor_det(M)

    @test detx(M') == conj(detx(M))


    MM = invx(M)
    @test M * MM == eye(Complex{BigInt}, 5)

    @test permanent(M') == permanent(M)'
end


@testset "Homogeneous" begin
    v = HVector(2, 3, 2)
    @test v[2] == 3 // 2
    @test dot(v, v) == 1
    w = HVector(20, 30, 20)
    @test v == w

    A = rand(Int, 3, 3) .% 20
    w = A * v
    @test w == HVector(A * v.data)
end

@testset "Permanent" begin
    A = ones(Int, 6, 6)
    @test permanent(A) == factorial(6)

    A = A - eye(Int, 6)
    @test permanent(A) == 265  # number of derangements of [6]

    A = triu(ones(5, 5))
    @test permanent(A) == 1
end

# @testset "Smith Normal Form" begin
#     A = rand(Mod{30}, 5, 5)
#     S, U, V = smith_normal_form(A)
#     @test S == U * A * V
#     @test isdiag(S)
#     @test all(LinearAlgebraX.divisible(S[i,i], S[i+1,i+1]) for i in 1:4)

#     B = rand(Mod{30}, 3, 3)
#     det_B = B[1,1]*B[2,2]*B[3,3] + B[1,2]*B[2,3]*B[3,1] + B[1,3]*B[2,1]*B[3,2] -
#             B[1,3]*B[2,2]*B[3,1] - B[1,1]*B[2,3]*B[3,2] - B[1,2]*B[2,1]*B[3,3]
#     smith_normal_form(B)
#     @test det_B == detx(B)
# end

# @testset "Triangular" begin
#     function det_lower_triangular(A::AbstractMatrix{T}) where {T<:AbstractMod}
#         L, col_ops = LinearAlgebraX.lower_triangular!(copy(A))
#         det_C = prod(detx(op) for op in col_ops)
#         # L = A * C
#         return prod(diag(L)) * inv(det_C)
#     end

#     A = rand(Mod{30}, 5, 5)
#     @test detx(A) == det_lower_triangular(A)
#     U, C = upper_triangular(A)
#     @test U == C * A
#     L, D = lower_triangular(A)
#     @test L == A * D
# end
