using Primes
export smith_normal_form

abstract type AbstractMatrixOperation{R} end
abstract type AbstractRowOperation{R} <: AbstractMatrixOperation{R} end
abstract type AbstractColumnOperation{R} <: AbstractMatrixOperation{R} end
struct RowSwap{R} <: AbstractRowOperation{R}
    i::Int
    j::Int
end
struct RowAddMult{R} <: AbstractRowOperation{R}
    i::Int
    j::Int
    c::R
end
struct RowScale{R} <: AbstractRowOperation{R}
    i::Int
    u::R
end
struct RowSmith{R} <: AbstractRowOperation{R}
    i::Int
    j::Int
    s::R
    u::R
    t::R
    v::R
end
struct ColumnSwap{R} <: AbstractColumnOperation{R}
    i::Int
    j::Int
end
struct ColumnAddMult{R} <: AbstractColumnOperation{R}
    i::Int
    j::Int
    c::R
end
struct ColumnScale{R} <: AbstractColumnOperation{R}
    i::Int
    u::R
end
struct ColumnSmith{R} <: AbstractColumnOperation{R}
    i::Int
    j::Int
    s::R
    u::R
    t::R
    v::R
end

detx(op::RowSwap{R}) where {R} = (op.i == op.j) ? one(R) : -one(R)
detx(::RowAddMult{R}) where {R} = one(R)
detx(op::RowScale{R}) where {R} = op.u
detx(op::RowSmith{R}) where {R} = op.s * op.v - op.u * op.t
detx(op::ColumnSwap{R}) where {R} = (op.i == op.j) ? one(R) : -one(R)
detx(::ColumnAddMult{R}) where {R} = one(R)
detx(op::ColumnScale{R}) where {R} = op.u
detx(op::ColumnSmith{R}) where {R} = op.s * op.v - op.u * op.t

(r::AbstractMatrixOperation{R})(M::AbstractMatrix{R}) where R = apply_matrix_operation!(M, r)

function apply_matrix_operation!(M::AbstractMatrix{R}, op::RowSwap{R}) where {R}
    M[op.i,:], M[op.j,:] = M[op.j,:], M[op.i,:]
    return M
end
function apply_matrix_operation!(M::AbstractMatrix{R}, op::RowAddMult{R}) where {R}
    M[op.j,:] += op.c * M[op.i,:]
    return M
end
function apply_matrix_operation!(M::AbstractMatrix{R}, op::RowScale{R}) where {R}
    M[op.i,:] *= op.u
    return M
end
function apply_matrix_operation!(M::AbstractMatrix{R}, op::RowSmith{R}) where {R}
    M[[op.i, op.j],:] = [op.s op.t; op.u op.v] * M[[op.i, op.j],:]
    return M
end
function apply_matrix_operation!(M::AbstractMatrix{R}, op::ColumnSwap{R}) where {R}
    M[:,op.i], M[:,op.j] = M[:,op.j], M[:,op.i]
    return M
end
function apply_matrix_operation!(M::AbstractMatrix{R}, op::ColumnAddMult{R}) where {R}
    M[:,op.j] += op.c * M[:,op.i]
    return M
end
function apply_matrix_operation!(M::AbstractMatrix{R}, op::ColumnScale{R}) where {R}
    M[:,op.i] *= op.u
    return M
end
function apply_matrix_operation!(M::AbstractMatrix{R}, op::ColumnSmith{R}) where {R}
    M[:,[op.i, op.j]] = M[:,[op.i, op.j]] * [op.s op.u; op.t op.v]
    return M
end

"""
smith_coeff(a::Mod{N, T}, b::Mod{N, T}[, prime_powers]) where {N, T}

Return a invertable matrix `M` such that `M*[a; b] = [g; 0]`, 
where `g` is a gcd of `a` and `b`.
"""
function smith_coeff(a::Mod{N, T}, b::Mod{N, T}, prime_powers::Union{Nothing, Vector{<:Integer}} = nothing) where {N, T}
    prime_powers = isnothing(prime_powers) ? [p^d for (p, d) in factor(N)] : prime_powers
    return smith_coeff(T, value(a), value(b), prime_powers, N)
end

function smith_coeff(T, a::Integer, b::Integer, prime_powers::Vector{<:Integer}, N::Integer)
    as = a .% prime_powers
    bs = b .% prime_powers
    zas = iszero.(as)
    zbs = iszero.(bs)
    za = Mod{N, T}(CRT(BigInt, zas, prime_powers))
    zb = Mod{N, T}(CRT(BigInt, zbs, prime_powers))
    zabs = zas .* zbs
    as = as + zas .* bs
    bs = bs + zbs .* as
    pas = gcd.(as, prime_powers)
    pbs = gcd.(bs, prime_powers)
    cas = as .÷ pas
    cbs = bs .÷ pbs
    pas = pas .* ((!).(zabs))
    pbs = pbs .* ((!).(zabs))
    cas = cas + zabs
    cbs = cbs + zabs
    @assert cas .* pas == as
    @assert cbs .* pbs == bs
    ca = Mod{N, T}(CRT(BigInt, cas, prime_powers))
    cb = Mod{N, T}(CRT(BigInt, cbs, prime_powers))
    αs = pbs .÷ (pas + zabs)
    βs = pas .÷ (pbs + zabs)
    γs = iszero.(αs)
    α = Mod{N, T}(CRT(BigInt, αs, prime_powers))
    β = Mod{N, T}(CRT(BigInt, βs, prime_powers))
    γ = Mod{N, T}(CRT(BigInt, γs, prime_powers))
    M_za = Mod{N, T}[1 za; 0 1]
    M_zb = Mod{N, T}[1 0; zb 1]
    M_inv_ca_cb = Mod{N, T}[inv(ca) 0; 0 inv(cb)]
    M_α = Mod{N, T}[1 0; -α 1]
    M_β = Mod{N, T}[1 -β; 0 1]
    M_plus_γ = Mod{N, T}[1 γ; 0 1]
    M_minus_γ = Mod{N, T}[1 -γ; 0 1]
    M = M_minus_γ * M_plus_γ * M_β * M_α * M_inv_ca_cb * M_zb * M_za
    res = Mod{N, T}[a; b]
    res = M * res
    @assert iszero(res[2])
    return M
end

"""
    divisible(a::Mod{N, T}, b::Mod{N, T}[, prime_powers]) where {N, T}

Check if a | b mod N, where N is the product of `prime_powers`.
"""
function divisible(a::Mod{N, T}, b::Mod{N, T}, prime_powers::Union{Nothing, Vector{<:Integer}} = nothing) where {N, T}
    prime_powers = isnothing(prime_powers) ? [p^d for (p, d) in factor(N)] : prime_powers
    _, r = mod_divrem(T, value(b), value(a), prime_powers, N)
    return iszero(r)
end

function Base.divrem(a::Mod{N, T}, b::Mod{N, T}, prime_powers::Union{Nothing, Vector{<:Integer}} = nothing) where {N, T}
    prime_powers = isnothing(prime_powers) ? [p^d for (p, d) in factor(N)] : prime_powers
    return mod_divrem(T, value(a), value(b), prime_powers, N)
end

function mod_divrem(T, a::Integer, b::Integer, prime_powers::Vector{<:Integer}, N::Integer)
    as = a .% prime_powers
    bs = b .% prime_powers
    pas = gcd.(as, prime_powers)
    pbs = gcd.(bs, prime_powers)

    # `a` is not divisible by `b`.
    !all(iszero, rem.(pas, pbs)) && return Mod{N, T}(0), Mod{N, T}(a)

    # `a` is divisible by `b`: `a = b * c`.
    zbs = iszero.(bs)
    αs = as .÷ pas 
    βs = bs .÷ pbs + zbs
    inv_βs = [gcdx(β, p)[2] for (β, p) in zip(βs, prime_powers)]
    cs = inv_βs .* αs .* (pas .÷ pbs)
    cs[zbs] .= 0
    # @assert all(rem.(ds, prime_powers) .== 0)
    c = Mod{N, T}(CRT(BigInt, cs, prime_powers))
    return c, Mod{N, T}(0)
end

"""
    smith_normal_form(M::Matrix{Mod{N, T}}[, prime_powers]) where {N, T}

Compute the Smith normal form of `M` over `Z_N`. It returns matrices `S`, `U`, and `V` 
such that `S = U * M * V` and `S` is in Smith normal form of `M`.
"""
function smith_normal_form(M::AbstractMatrix{Mod{N, T}}, prime_powers::Union{Nothing, Vector{<:Integer}} = nothing) where {N, T}
    m, n = size(M)
    S, row_ops, col_ops = smith_normal_form!(copy(M), prime_powers)
    U = eye(Mod{N, T}, m)
    V = eye(Mod{N, T}, n)
    for op in row_ops
        op(U)
    end
    for op in col_ops
        op(V)
    end
    return S, U, V
end

function smith_normal_form!(M::AbstractMatrix{Mod{N, T}}, prime_powers::Union{Nothing, Vector{<:Integer}} = nothing) where {N, T}
    prime_powers = isnothing(prime_powers) ? [p^d for (p, d) in factor(N)] : prime_powers
    R = Mod{N, T}
    m, n = size(M)
    row_ops = AbstractMatrixOperation{R}[]
    col_ops = AbstractMatrixOperation{R}[]
    for i = 1:min(m, n)
        smith_elimination!(M, row_ops, col_ops, i, prime_powers)
    end
    diag_normalized = false
    while diag_normalized
        diag_normalized = true
        for i = 2:min(m, n)
            if !divisible(M[i-1,i-1], M[i,i])
                diag_normalized = false
                op = ColumnAddMult(i, i-1, one(R))
                push!(col_ops, op)
                op(M)
                smith_elimination!(M, row_ops, col_ops, i-1, prime_powers)
                break
            end
        end
    end
    return M, row_ops, col_ops
end

function smith_elimination!(M::AbstractMatrix{Mod{N, T}}, row_ops, col_ops, i, prime_powers::Vector{<:Integer}) where {N, T}
    R = Mod{N, T}
    m, n = size(M)
    @assert 1 <= i <= min(m, n)
    finished = false
    while !finished
        finished = true
        
        # Eliminate the i-th column.
        if iszero(M[i,i])
            for j = i+1:m
                if !iszero(M[j,i])
                    op = RowSwap{R}(i, j)
                    push!(row_ops, op)
                    op(M)
                    break
                end
            end
        end
        for j = i+1:m
            c, r = divrem(M[j,i], M[i,i])
            if iszero(r)
                op = RowAddMult{R}(i, j, -c)
                push!(row_ops, op)
                op(M)
            else
                op = RowSmith{R}(i, j, smith_coeff(M[i,i], M[j,i], prime_powers)[:]...)
                push!(row_ops, op)
                op(M)
                finished = false
            end
        end

        # Eliminate the i-th row.
        if iszero(M[i,i])
            for j = i+1:n
                if !iszero(M[i,j])
                    op = ColumnSwap{R}(i, j)
                    push!(col_ops, op)
                    op(M)
                    break
                end
            end
        end
        for j = i+1:n
            c, r = divrem(M[i,j], M[i,i])
            if iszero(r)
                op = ColumnAddMult{R}(i, j, -c)
                push!(col_ops, op)
                op(M)
            else
                op = ColumnSmith{R}(i, j, smith_coeff(M[i,i], M[i,j], prime_powers)[:]...)
                push!(col_ops, op)
                op(M)
                finished = false
            end
        end
    end
    return M, row_ops, col_ops
end
