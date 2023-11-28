export eye

"""
    eye(T, n)::Matrix

`eye(T,n)` returns an `n`-by-`n` identity matrix whose entries
are numbers of type `T`.

`eye(n)` is the same as `eye(Float64,n)`.
"""
function eye(T, n)::Matrix
    return Matrix{T}(I, n, n)
end

eye(n) = eye(Float64, n)
