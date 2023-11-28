export detx

cofactor_warning = "Using cofactor expansion to calculate determinant; may be very slow."

"""
    detx(A::AbstractMatrix{T}) where {T<:IntegerX}

`detx(A::AbstractMatrix{T})` is an exact determinant of the matrix.
Here `T` can be any kind of integer, rational, or `Mod`.
"""
function detx(A::AbstractMatrix{T}) where {T<:Union{IntegerX,RationalX}}
    A = big.(A)
    TT = typeof(A[1, 1])
    return TT(det(A // 1))
end

function detx(A::AbstractMatrix{T})::T where {T<:Mod}
    try
        return det(A)
    catch
        @info cofactor_warning
        return cofactor_det(A)
    end

end

# Everything else
function detx(A::AbstractMatrix{T})::T where {T}
    B = collect(A)
    try
        return detx!(B)
    catch
        @info cofactor_warning
        return cofactor_det(B)  # if all else fails!
    end
end


# detx! is the work behind det when we can't use Julia's det.
# it can modify the matrix so this is not exposed to the general public.

function detx!(A::AbstractMatrix{T})::T where {T}
    r = LinearAlgebra.checksquare(A)

    if r == 0
        return one(T)
    end

    if r == 1
        return A[1, 1]
    end

    if r == 2
        return A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
    end

    col = A[:, 1]  # first column
    if all(col .== 0)
        return zero(T)
    end

    k = findfirst(col .!= 0)  # we pivot here

    sign_factor = one(Int8)
    if k > 1
        sign_factor = -sign_factor
        row_swap!(A, 1, k)
    end

    factor = A[1, 1] * sign_factor  # multiply by this at the end

    row_scale!(A, 1, 1 // A[1, 1])

    for i = 2:r
        row_add_mult!(A, 1, -A[i, 1], i)
    end


    return detx!(A[2:end, 2:end]) * factor

end

function detx!(A::AbstractMatrix{<:Mod})
    U, row_ops = upper_triangular!(A)
    det_C = prod(detx(op) for op in row_ops)

    # U = C * A
    return inv(det_C) * prod(diag(U))
end