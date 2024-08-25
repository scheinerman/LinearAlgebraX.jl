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
