"""
`hilbert(n)` creates an `n`-by-`n` Hilbert matrix.

## Example
```
julia> hilbert(4)
4×4 Array{Rational{Int64},2}:
 1//1  1//2  1//3  1//4
 1//2  1//3  1//4  1//5
 1//3  1//4  1//5  1//6
 1//4  1//5  1//6  1//7
```
"""
function hilbert(n::Int)
    A = zeros(Rational{Int}, n, n)
    for i = 1:n
        for j = 1:n
            A[i, j] = 1 // (i + j - 1)
        end
    end
    return A
end
