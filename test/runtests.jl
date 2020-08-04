using Test
using LinearAlgebra,LinearAlgebraX, Mods

A = big.(rand(Int,10,10) .% 100)

include("hilbert.jl")

H = hilbert(12)
@test rankx(H) == 12

@test detx(A) == det(Rational.(A))
@test detx(A) == detx(A')
@test A * invx(A) == eye(BigInt, 10)

A = ones(Int,3,5)
N = nullspacex(A)
@test iszero(A*N)


using SimplePolynomials

A = triu(ones(Int,5,5))
p = char_poly(A)
x = getx()
@test p == (x-1)^5

A = ones(Int,5,5)
p = char_poly(A)
@test p == x^4 * (x-5)

T = Mod{17}
A = rand(T, 6, 6)
p = char_poly(A)
@test detx(A) == p(0)


A = rand(Int,20,20) .% 100
@test rankx(A) <= 20
