export permanent

_enlarge(x::T) where {T<:Union{IntegerX,RationalX}} = big(x)
_enlarge(x) = x

"""
    permanent(A)

Compute the permanent of the square matrix `A`.
"""
function permanent(A::AbstractMatrix{T}) where {T}     # code by Daniel Scheinerman
    A = _enlarge.(A)
    m, n = size(A)
    s = _enlarge(zero(T)) #to accumulate the permanent
    x = UInt64(0) #gray code state
    v = zeros(Int64, n) #current subset sum
    for i = 1:2^n-1
        t = trailing_zeros(i) #next bit flip
        b = (x >> t) & 1 #if 1 then add row, else subtract
        x = xor(x, 1 << t)
        if b == 0
            v += A[t+1, :]
        else
            v -= A[t+1, :]
        end
        s += (-1)^(i & 1) * prod(v)
    end
    s *= (-1)^n
    return s
end


# Old code available as LinearAlgebraX.old_permanent

function old_permanent(A::AbstractMatrix{T}) where {T<:Union{IntegerX,RationalX}}
    A = big.(A)
    return perm_work(A)
end

function old_permanent(A::AbstractMatrix{T})::T where {T}
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
