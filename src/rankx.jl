export rankx

"""
`rankx(A)` computes the rank of the exact matrix `A`
"""
function rankx(A::AbstractMatrix{T})::Int where {T}
    r, c = size(A)

    AA = deepcopy(A)
    try
        AA //= 1  # convert to rational if need be
        try
            rrefx!(AA)
        catch
            AA = big.(AA)
            rrefx!(AA)
        end
    catch
        rrefx!(AA)
    end

    count = 0
    for i = 1:r
        if any(AA[i, :] .!= 0)
            count += 1
        end
    end
    return count
end
