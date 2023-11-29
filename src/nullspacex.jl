export nullspacex

"""
`nullspacex(A)` returns an exact basis for the nullspace of the matrix `A`
"""
nullspacex(A::AbstractMatrix{T}) where T = _nullspacex(A)


function _nullspacex(A::AbstractMatrix{T}) where {T}
    r, c = size(A)
    B = rrefx(A)

    leads = Int[]
    # in each row, find first 1
    for i = 1:r
        row = B[i, :]
        if all(row .== 0)
            continue
        end
        k = findfirst(row .!= 0)
        append!(leads, k)
    end

    frees = setdiff(collect(1:c), leads)


    result = Matrix{T}(undef, c, 0)
    for f in frees
        v = zeros(T, c)
        v[f] = T(1)
        for i = 1:length(leads)
            l = leads[i]  # (i,l) is a leading 1
            if l < f
                v[l] = -B[i, f]
            end
        end
        result = [result v]
    end
    return result
end


function nullspacex(A::AbstractMatrix{T}) where {T<:IntegerX}
    return nullspacex(big.(A) // 1)
end

function nullspacex(A::AbstractMatrix{T}) where T<:AbstractAlgebraicFunction
    A = SimpleRationalFunction.(A)
    return _nullspacex(A)
end