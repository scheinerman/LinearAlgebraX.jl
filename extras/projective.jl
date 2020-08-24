using Mods, Primes, SimpleGraphs, LinearAlgebraX

"""
`generate_points(p::Int)` generates the homogeneous coordinates of the 
points in a finite projective plane of prime order `p`.
"""
function generate_points(p::Integer)
    if ~isprime(p)
        error("Modulus $p is not prime")
    end

    n = p * p + p + 1
    A = Array{HVector}(undef, n)
    count = 1

    # case c=1
    for a = 0:p-1
        for b = 0:p-1
            x = Mod{p}.([a, b, 1])
            A[count] = HVector(x)
            count += 1
        end
    end

    # case c=0, b=1
    for a = 0:p-1
        x = Mod{p}.([a, 1, 0])
        A[count] = HVector(x)
        count += 1
    end

    A[count] = HVector(Mod{p}(1), 0, 0)

    return A
end

"""
`incidence_matrix(p::Int)` creates the point/line incidence matrix 
of a finite projective plane of prime order `p`.
"""
function incidence_matrix(p::Int)
    pts = generate_points(p)
    n = p*p + p + 1
    A = zeros(Int,n,n)
    for u = 1:n
        P = pts[u]
        for v = 1:n
            Q = pts[v]
            if dot(P,Q)==1
                A[u,v] = 1
            end
        end
    end
    return A
end

"""
`incidence_graph(p::Int)` creates the bipartite point-line 
incidence graph of a finite projective plane of prime order `p`.
"""
function incidence_graph(p::Int)
    M = incidence_matrix(p)
    r,c = size(M)
    A = [ zeros(r,c)  M; copy(M')  zeros(r,c) ]
    G = SimpleGraph(A)
end 