export upper_triangular, lower_triangular

function upper_triangular(M::AbstractMatrix{Mod{N, T}}, prime_powers::Union{Nothing, Vector{<:Integer}} = nothing) where {N, T}
    U = copy(M)
    U, row_ops = upper_triangular!(U, prime_powers)
    C = eye(Mod{N, T}, size(U, 1))
    for op in row_ops
        apply_matrix_operation!(C, op)
    end
    return U, C
end

function upper_triangular!(M::AbstractMatrix{Mod{N, T}}, prime_powers::Union{Nothing, Vector{<:Integer}} = nothing) where {N, T}
    prime_powers = isnothing(prime_powers) ? [p^d for (p, d) in factor(N)] : prime_powers
    R = Mod{N, T}
    ops = AbstractMatrixOperation{R}[]
    m, n = size(M)
    for i = 1:min(m, n)
        # Find a non-zero element in the i-th column.
        if iszero(M[i,i])
            for j = i+1:m
                if !iszero(M[j,i])
                    op = RowSWAP{R}(i, j)
                    push!(ops, op)
                    apply_matrix_operation!(M, op)
                    break
                end
            end
        end
        # Make U[i,j] = 0 for j = i+1:n.
        for j = i+1:m
            if !iszero(M[j,i])
                c, r = divrem(M[j,i], M[i,i], prime_powers)
                op = if iszero(r)
                    RowAddTo{R}(i, j, -c)
                else
                    RowSmith{R}(i, j, smith_coeff(M[i,i], M[j,i], prime_powers)[:]...)
                end
                @show op
                push!(ops, op)
                apply_matrix_operation!(M, op)
            end
        end
    end
    return M, ops
end

function lower_triangular(M::AbstractMatrix{Mod{N, T}}, prime_powers::Union{Nothing, Vector{<:Integer}} = nothing) where {N, T}
    L = copy(M)
    L, col_ops = lower_triangular!(L, prime_powers)
    C = eye(Mod{N, T}, size(L, 2))
    for op in col_ops
        apply_matrix_operation!(C, op)
    end
    return L, C
end
function lower_triangular!(M::AbstractMatrix{Mod{N, T}}, prime_powers::Union{Nothing, Vector{<:Integer}} = nothing) where {N, T}
    prime_powers = isnothing(prime_powers) ? [p^d for (p, d) in factor(N)] : prime_powers
    R = Mod{N, T}
    ops = AbstractMatrixOperation{R}[]
    m, n = size(M)
    for i = 1:min(m, n)
        # Find a non-zero element in the -th row.
        if iszero(M[i,i])
            for j = i+1:n
                if !iszero(M[i,j])
                    op = ColumnSWAP{R}(i, j)
                    push!(ops, op)
                    apply_matrix_operation!(M, op)
                    break
                end
            end
        end
        # Make U[i,j] = 0 for j = i+1:n.
        for j = i+1:n
            if !iszero(M[i,j])
                c, r = divrem(M[i,j], M[i,i], prime_powers)
                op = if iszero(r)
                    ColumnAddTo{R}(i, j, -c)
                else
                    ColumnSmith{R}(i, j, smith_coeff(M[i,i], M[i,j], prime_powers)[:]...)
                end
                push!(ops, op)
                apply_matrix_operation!(M, op)
            end
        end
    end
    return M, ops
end
