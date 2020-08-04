export rankx

"""
`rankx(A)` computes the rank of the exact matrix `A`
"""
function rankx(A::AbstractArray{T,2})::Int where T
    r,c = size(A)

    AA = A//1
    try
        rrefx!(AA)
    catch
        AA = big.(AA)
        rrefx!(AA)
    end


    count = 0
    for i=1:r
        if !iszero(AA[i,:])
            count += 1
        end
    end
    return count
end


"""
`nullspacex(A)` returns an exact basis for the matrix `A`
"""
function nullspacex(A::AbstractArray{T,2}) where T
    r,c = size(A)
    B = rrefx(A)

    leads = Int[]
    # in each row, find first 1
    for i=1:r
        row = B[i,:]
        k = findfirst(row .!= 0)
        if k !== nothing
            append!(leads,k)
        end
    end

    frees = setdiff(collect(1:c), leads)

    result = Matrix{T}(c,0)
    for j in frees
        v = zeros(T,c)
        v[j] = 1
        for i in leads
            if i<j
                v[i] = -B[i,j]
            end
        end
        result = [result v]
    end

    return result
end
