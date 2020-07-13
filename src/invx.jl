# matrix (and other) inversion

export invx

"""
`invx(A)` for a matrix `A` gives an exact matrix inverse.
"""
function invx(A::AbstractArray{T,2}) where T
    r,c = size(A)
    @assert r==c "Matrix must be square"

    if r==0
        return A
    end

    B = eye(T,r)

    X = [A B]
    rrefx!(X)

    d = [X[i,i] for i=1:r]
    if all(d.==1)
        return X[:,r+1:end]
    end
    @error "Matrix is not intertible"
end
