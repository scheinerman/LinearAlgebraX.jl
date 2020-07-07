
# IntegerX is any sort of real or Gaussian integer
IntegerX = Union{S,Complex{S}} where S<:Integer

function _recip(x::T) where T <: IntegerX
    return 1//x
end
_recip(x) = inv(x)


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
    end


    col = A[:,1]  # first column
    if all(col .== 0)
        return 0
    end

    k = findfirst(col .!== 0)  # we pivot here
    b = _recip(A[k,1])
    # row = b * A[k,:]

    row = [b * A[k,j] for j=1:c]
    row = reshape(row,1,r)

    println("row = $row")


    for i=1:r
        if i != k  # leave row k alone
            for j=1:c
                A[i,j] -= b*A[i,1]*A[i,j]
            end

            #A[i,:] = A[i,:] -  A[i,1].*row
        end
    end

    idx = [collect(1:k-1); collect(k+1:r)]

    B = A[idx,2:end]

    return (-1)^(k-1) * A[k,1] * detx(B)
end
