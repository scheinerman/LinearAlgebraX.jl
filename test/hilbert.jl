function hilbert(n::Int)
    A = zeros(Rational{Int}, n, n)
    for i = 1:n
        for j = 1:n
            A[i, j] = 1 // (i + j - 1)
        end
    end
    return A
end
