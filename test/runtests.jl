using Test
using LinearAlgebra,LinearAlgebraX

A = big.(rand(Int,10,10) .% 100)

include("hilbert.jl")

H = hilbert(12)
@test rankx(H) == 12

@test detx(A) == det(Rational.(A))
@test detx(A) == detx(A')
@test A * invx(A) == eye(BigInt, 10)

A = ones(Int,3,5)
N = nullspacex(A)
@test all(0 .== A*N)
