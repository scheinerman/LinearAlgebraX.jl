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

* `detx` -- exact determinant
* `nullspacex` -- exact nullspace
* `rankx` -- exact rankx
* `invx` -- exact inverse

#### Example

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
