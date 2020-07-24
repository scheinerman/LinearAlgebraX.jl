export detx

"""
`detx(A::Matrix{T})` is an exact determinant of the matrix.
Here `T` can be any kind of integer, rational, `Mod`, or `GF2`.

I hope to expand this to `Polynomial`s.
"""
function detx(A::AbstractArray{T,2}) where T<:IntegerX
    # @info "Using IntegerX detx{$T}"
    try
        return T(det(A//1))
    catch
        A = big.(A)
        return detx(A)
    end
end

function detx(A::AbstractArray{T,2}) where T<:RationalX
    # @info "Using RationalX detx{$T}"
    try
        return det(A)
    catch
        return det(big.(A))
    end
end

function detx(A::AbstractArray{T,2}) where T
    # @info "Using detx! on matrix of type $T"
    r,c = size(A)
    B = Matrix{Any}(undef,r,c)
    for i=1:r
        for j=1:c
            B[i,j] = A[i,j]
        end
    end
    # B = copy(A)
    return detx!(B)
end


# detx! is the work behind det when we can't use Julia's det.
# it can modify the matrix so this is not exposed to the general public.

function detx!(A::AbstractArray{T,2}) where T
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

    col = A[:,1]  # first column
    if all(col .== 0)
        return 0
    end

    k = findfirst(col .!= 0)  # we pivot here

    sign_factor = 1
    if k>1
        sign_factor = -sign_factor
        row_swap!(A,1,k)
    end

    factor = A[1,1] * sign_factor  # multiply by this at the end

    row_scale!(A,1,1//A[1,1])

    for i=2:r
        row_add_mult!(A,1,-A[i,1],i)
    end


    return detx!(A[2:end,2:end]) * factor

end
