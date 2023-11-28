export permanent
"""
    permanent(A)

Compute the permanent of the square matrix `A`. Works best if `A`
is relatively sparse. For a dense matrix, this will be slow.
"""
function permanent(A::AbstractMatrix{T}) where {T<:Union{IntegerX,RationalX}}
    A = big.(A)
    return perm_work(A)
end

function permanent(A::AbstractMatrix{T})::T where {T}
    perm_work(A)
end

function perm_work(A::AbstractMatrix{T}) where {T}
    r, c = size(A)
    if r != c
        error("Matrix must be square")
    end
    if r == 0
        return zero(T)
    end

    allow = [findall((x) -> x != 0, A[i, :]) for i = 1:r]
    PG = PermGen(allow)

    try
        first(PG)  # make sure PG is not empty
        return sum(prod(A[i, p[i]] for i = 1:r) for p in PG)
    catch
        return zero(T)
    end
end
