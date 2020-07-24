# LinearAlgebraX



[![Build Status](https://travis-ci.org/scheinerman/LinearAlgebraX.jl.svg?branch=master)](https://travis-ci.org/scheinerman/LinearAlgebraX.jl)

[![Coverage Status](https://coveralls.io/repos/scheinerman/LinearAlgebraX.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/scheinerman/LinearAlgebraX.jl?branch=master)

[![codecov.io](http://codecov.io/github/scheinerman/LinearAlgebraX.jl/coverage.svg?branch=master)](http://codecov.io/github/scheinerman/LinearAlgebraX.jl?branch=master)


This module implements basic linear algebra methods for matrices
with exact entries (e.g., `Rational{Int}` values).  The function names
typically match the standard ones in Julia but with an `x` (for "exact")
appended.



## Functions

These functions in this module end with the letter `x`
and have the same definitions as their counterparts that do not have an `x`.
For exact types (such as `Int`s) these functions give exact results.

* `detx` -- exact determinant (via row reduced echelon form)
* `cofactor_det`-- slower exact determinant (via cofactor expansion)
* `nullspacex` -- exact nullspace
* `rankx` -- exact rankx
* `invx` -- exact inverse
* `rrefx` -- row reduced echelon form
* `eye` -- lovingly restored
* `char_poly` -- characteristic polynomial

## Examples

#### Determinant
```
julia> A = ones(Int,10,10)+eye(Int,10);

julia> det(A)
11.000000000000004

julia> detx(A)
11

julia> A = rand(Int,20,20) .% 20;

julia> det(A)
3.3905496651565455e29

julia> detx(A)
339054966515654744413389494504
```

#### Nullspace

```
julia> A = reshape(collect(1:12),3,4)
3×4 Array{Int64,2}:
 1  4  7  10
 2  5  8  11
 3  6  9  12

julia> nullspacex(A)
4×2 Array{Rational{BigInt},2}:
  1//1   2//1
 -2//1  -3//1
  1//1   0//1
  0//1   1//1

julia> nullspace(A)
4×2 Array{Float64,2}:
 -0.475185  -0.272395
  0.430549   0.717376
  0.564458  -0.617566
 -0.519821   0.172585
```

#### Rank

Consider the 12-by-12 Hibert matrix, `H`.
```
12×12 Array{Rational{Int64},2}:
 1//1   1//2   1//3   1//4   1//5   1//6   1//7   1//8   1//9   1//10  1//11  1//12
 1//2   1//3   1//4   1//5   1//6   1//7   1//8   1//9   1//10  1//11  1//12  1//13
 1//3   1//4   1//5   1//6   1//7   1//8   1//9   1//10  1//11  1//12  1//13  1//14
 1//4   1//5   1//6   1//7   1//8   1//9   1//10  1//11  1//12  1//13  1//14  1//15
 1//5   1//6   1//7   1//8   1//9   1//10  1//11  1//12  1//13  1//14  1//15  1//16
 1//6   1//7   1//8   1//9   1//10  1//11  1//12  1//13  1//14  1//15  1//16  1//17
 1//7   1//8   1//9   1//10  1//11  1//12  1//13  1//14  1//15  1//16  1//17  1//18
 1//8   1//9   1//10  1//11  1//12  1//13  1//14  1//15  1//16  1//17  1//18  1//19
 1//9   1//10  1//11  1//12  1//13  1//14  1//15  1//16  1//17  1//18  1//19  1//20
 1//10  1//11  1//12  1//13  1//14  1//15  1//16  1//17  1//18  1//19  1//20  1//21
 1//11  1//12  1//13  1//14  1//15  1//16  1//17  1//18  1//19  1//20  1//21  1//22
 1//12  1//13  1//14  1//15  1//16  1//17  1//18  1//19  1//20  1//21  1//22  1//23
```
We compare the results of `rank` (from the `LinearAlgebra` module) and
`rankx` (in this module):
```
julia> rank(H)
11

julia> rankx(H)
12
```

#### Inverse

```
julia> using Mods

julia> A = rand(Mod{11},5,5)
5×5 Array{Mod{11},2}:
 Mod{11}(2)   Mod{11}(4)  Mod{11}(4)  Mod{11}(0)   Mod{11}(2)
 Mod{11}(9)   Mod{11}(4)  Mod{11}(5)  Mod{11}(1)  Mod{11}(10)
 Mod{11}(3)   Mod{11}(4)  Mod{11}(5)  Mod{11}(6)   Mod{11}(0)
 Mod{11}(5)  Mod{11}(10)  Mod{11}(4)  Mod{11}(5)   Mod{11}(4)
 Mod{11}(9)  Mod{11}(10)  Mod{11}(7)  Mod{11}(8)   Mod{11}(9)

julia> B = invx(A)
5×5 Array{Mod{11},2}:
 Mod{11}(4)  Mod{11}(5)  Mod{11}(0)   Mod{11}(6)   Mod{11}(8)
 Mod{11}(7)  Mod{11}(4)  Mod{11}(9)  Mod{11}(10)   Mod{11}(3)
 Mod{11}(6)  Mod{11}(0)  Mod{11}(2)   Mod{11}(5)   Mod{11}(5)
 Mod{11}(3)  Mod{11}(4)  Mod{11}(9)  Mod{11}(10)  Mod{11}(10)
 Mod{11}(9)  Mod{11}(9)  Mod{11}(0)   Mod{11}(8)   Mod{11}(9)

julia> A*B
5×5 Array{Mod{11},2}:
 Mod{11}(1)  Mod{11}(0)  Mod{11}(0)  Mod{11}(0)  Mod{11}(0)
 Mod{11}(0)  Mod{11}(1)  Mod{11}(0)  Mod{11}(0)  Mod{11}(0)
 Mod{11}(0)  Mod{11}(0)  Mod{11}(1)  Mod{11}(0)  Mod{11}(0)
 Mod{11}(0)  Mod{11}(0)  Mod{11}(0)  Mod{11}(1)  Mod{11}(0)
 Mod{11}(0)  Mod{11}(0)  Mod{11}(0)  Mod{11}(0)  Mod{11}(1)
 ```

 #### Characteristic polynomial

```
julia> using SimplePolynomials, LinearAlgebra

julia> x = getx()
x

julia> A = triu(ones(Int,5,5))
5×5 Array{Int64,2}:
 1  1  1  1  1
 0  1  1  1  1
 0  0  1  1  1
 0  0  0  1  1
 0  0  0  0  1

julia> char_poly(A)
-1 + 5*x - 10*x^2 + 10*x^3 - 5*x^4 + x^5

julia> ans == (x-1)^5
true

julia> using Mods

julia> A = rand(Mod{17},4,4)
4×4 Array{Mod{17},2}:
 Mod{17}(16)  Mod{17}(10)   Mod{17}(9)  Mod{17}(12)
 Mod{17}(15)   Mod{17}(1)   Mod{17}(1)   Mod{17}(6)
  Mod{17}(3)   Mod{17}(2)   Mod{17}(5)  Mod{17}(11)
  Mod{17}(5)  Mod{17}(15)  Mod{17}(15)   Mod{17}(7)

julia> char_poly(A)
Mod{17}(1) + Mod{17}(1)*x + Mod{17}(16)*x^2 + Mod{17}(5)*x^3 + Mod{17}(1)*x^4

julia> detx(A)
Mod{17}(1)
```

 #### Row reduced echelon form

 ```
 julia> A = rand(Int,4,6) .% 10
4×6 Array{Int64,2}:
 6   8  0  -6  -5   4
 0  -5  2   0  -3  -4
 0  -4  2  -8   7  -8
 1  -3  7   2  -6   2

julia> c = A[:,1] + A[:,2] - A[:,3]
4-element Array{Int64,1}:
 14
 -7
 -6
 -9

julia> A = [c A]
4×7 Array{Int64,2}:
 14  6   8  0  -6  -5   4
 -7  0  -5  2   0  -3  -4
 -6  0  -4  2  -8   7  -8
 -9  1  -3  7   2  -6   2

julia> rrefx(A)
4×7 Array{Rational{Int64},2}:
 1//1  0//1  0//1  -1//1  0//1   -23//130  -36//65
 0//1  1//1  0//1   1//1  0//1  -883//325  158//325
 0//1  0//1  1//1   1//1  0//1   551//650  512//325
 0//1  0//1  0//1   0//1  1//1  -379//325  204//325
 ```

 ## To do

 Still having some issues with integer overflow.
