using Test
using LinearAlgebra,LinearAlgebraX

A = rand(Int,6,6) .% 100

@test detx(A) == det(Rational.(A))


@test true
