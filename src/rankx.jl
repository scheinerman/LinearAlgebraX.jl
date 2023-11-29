export rankx

"""
    rankx(A::AbstractMatrix{T})::Int where {T}

`rankx(A)` computes the rank of the exact matrix `A`
"""
function rankx(A::AbstractMatrix{T})::Int where T
    r,c = size(A)
    return c-nullityx(A)
end
