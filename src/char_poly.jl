export char_poly

function char_poly(A::Matrix{T}) where T<:TypeX
    r,c = size(A)
    @assert r==c "Matrix must be square"

    xI = zeros(SimplePolynomial,r,r)
    x = getx()
    for i=1:r
        xI[i,i] = x
    end

    f = detx(xI-A)
    p = integerize(numerator(f))
end
