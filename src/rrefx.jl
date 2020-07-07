function rrefx(A::AbstractArray{T,2}) where T
    AA = copy(A)
    return rrefx!(AA)
end


function rrefx!(A::AbstractArray{T,2}) where T
    r,c = size(A)
    if c==0
        return A
    end

    # if the first column is all 0s, delete and move on
    col = A[:,1]
    if all(col .== 0)
        B = rrefx!(A[:,2:end])
        return [col B]
    end

    # get the location of the 1st nonzero entry in column 1
    k = findfirst(col .!= 0)

    # swap rows 1 and k, and make 1st entry a 1

    a = A[k,1]
    row = [ A[k,j]/a for j=1:c ]

    # (1/a) .* A[k,:]
    A[k,:] = A[1,:]
    A[1,:] = row

    # if there's only one row, we're done
    if r==1
        return A
    end

    # use the [1,1] entry to clean out all the rows
    for i=2:r
        b = A[i,1]
        for j=1:c
            A[i,j] = A[i,j] - b*row[j]
        end
    end

    # if only one column, we're done
    if c==1
        return A
    end

    col = A[:,1]
    row = A[1,:]


end
