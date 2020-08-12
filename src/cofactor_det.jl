export cofactor_det

"""
`cofactor(A::AbstractMatrix{T})` computes the determinant of `A` using
cofactor expansion (which can be slow).
The return type of this method is a number of type `T`.
The entires in `A` can be polynomials and that won't work with
Julia's `det`.
"""
function cofactor_det(A::AbstractMatrix{T}) where T
    r,c = size(A)
    @assert r==c "Matrix must be square"
    @assert r>0 "Matrix cannot be 0-by-0"
    if r==1
        return A[1,1]
    end

    total = 0 * A[1,1]
    sign = -1
    for i=1:r
        a = A[i,1]
        sign = -sign
        if a==0
            continue
        end
        ii = vcat(1:i-1, i+1:r)
        jj = 2:r
        B = A[ii,jj]
        total += sign * a * cofactor_det(B)
    end

    return total
end
