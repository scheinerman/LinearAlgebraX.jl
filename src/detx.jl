

# IntegerX is any sort of real or Gaussian integer
IntegerX = Union{S,Complex{S}} where S<:Integer

function _recip(x::T) where T <: IntegerX
    return 1//x
end
_recip(x) = inv(x)

"""
`detx(A::Matrix{T})` is an exact determinant of the matrix.
Here `T` can be any kind of integer, rational, `Mod`, or `GF2`.

I hope to expand this to `Polynomial`s.
"""
function detx(A::Matrix{T}) where T
    r,c = size(A)
    @assert r==c "Matrix must be square"

    if r==0
        return 1
    end

    if r==1
        return A[1,1]
    end

    if r==2
        return A[1,1]*A[2,2] - A[1,2]*A[2,1]
    end

    if T <: IntegerX
        A = _recip(one(T)) .* A  # make it rational
    else
        A = copy(A)
    end


    col = A[:,1]  # first column
    if all(col .== 0)
        return 0
    end

    k = findfirst(col .!= 0)  # we pivot here

    b = _recip(A[k,1])
    row = [b * A[k,j] for j=1:c]


    for i=1:r
        if i != k  # leave row k unchanged
            a = A[i,1]
            for j=1:c
                A[i,j] = A[i,j] - a * row[j]
            end
        end
    end

    idx = [collect(1:k-1); collect(k+1:r)]
    B = A[idx,2:end]

    d =  (-1)^(k-1) * A[k,1] * detx(B)

    if T <: IntegerX
        return T(d)
    end

    return d
end
