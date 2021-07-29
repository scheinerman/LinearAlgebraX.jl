# Elementary row operations

export row_swap!, row_scale!, row_add_mult!


"""
`row_swap!(A,i,j)` swaps rows `i` and `j` in the matrix `A`.
"""
function row_swap!(A::AbstractMatrix{T}, i::Int, j::Int) where {T}
    r, c = size(A)
    @assert (1 <= i <= r) && (1 <= j <= r) "Row index out of bounds"
    if i == j
        return
    end

    save_row = [A[i, k] for k = 1:c]

    for k = 1:c
        @inbounds A[i, k] = A[j, k]
    end

    for k = 1:c
        @inbounds A[j, k] = save_row[k]
    end
end

"""
`row_scale!(A,i,s)` multiplies all entries in row `i`
of `A` by `s`.
"""
function row_scale!(A::AbstractMatrix{T}, i::Int, s) where {T}
    r, c = size(A)
    @assert 1 <= i <= r "Row index out of bounds"

    for k = 1:c
        @inbounds A[i, k] = s * A[i, k]
    end
end


"""
`row_add_mult!(A,i,s,j)` adds `s` times row `i` to row `j`
in the matrix `A`.
"""
function row_add_mult!(A::AbstractMatrix{T}, i::Int, s, j::Int) where {T}
    r, c = size(A)
    @assert (1 <= i <= r) && (1 <= j <= r) "Row index out of bounds"

    save_row = [s * A[i, k] for k = 1:c]
    for k = 1:c
        @inbounds A[j, k] = A[j, k] + save_row[k]
    end
end
