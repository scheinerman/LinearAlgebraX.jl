export char_poly

"""
    char_poly(A::AbstractMatrix{T}) where {T<:TypeX}

`char_poly(A)` returns the characteristic polynomial of
the exact matrix `A`.
"""
function char_poly(A::AbstractMatrix{T}) where {T<:TypeX}
    r, c = size(A)
    @assert r == c "Matrix must be square"

    xI = zeros(SimplePolynomial, r, r)
    x = getx(T)
    for i = 1:r
        xI[i, i] = x
    end

    f = detx(xI - A)
    p = integerize(numerator(f))
end
