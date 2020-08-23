using Test
using LinearAlgebra, LinearAlgebraX, Mods


@testset "Basics" begin
    A = big.(rand(Int, 10, 10) .% 100)

    include("hilbert.jl")

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


@testset "Homogeneous" begin 
    v = HVector(2,3,2)
    @test v[2] == 3//2
    @test dot(v,v)==1
    w = HVector(20,30,20)
    @test v==w 

    A = rand(Int,3,3) .% 20
    w = A*v
    @test w == HVector(A*v.data)
end