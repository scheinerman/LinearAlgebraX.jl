# Extras for `LinearAlgebraX`

## Projective planes

The file `projective.jl` contains code for creating structures related 
to finite projective planes of prime order.

* `generate_points(p)` generates a list of the homogeneous coordinates of the 
points (lines) of a finite projective plane of prime order `p`.
* `incidence_matrix(p)` generates the point-line incidence matrix of a 
projective plane of order `p`.
* `incidence_graph(p)` generates the bipartite point-line incidence graph 
of a finite projective plane of order `p`.
* `incidence_hypergraph(p)` generates a hypergraph whose vertices are the 
points of an order-`p` projective plane and whose hyperedges are the lines 
of that plane. 
```julia
julia> include("extras/projective.jl")
incidence_graph

julia> A = incidence_matrix(3)
13×13 Array{Int64,2}:
 0  0  0  0  0  0  0  0  0  1  1  1  1
 0  0  1  0  0  1  0  0  1  0  0  0  1
 0  1  0  0  1  0  0  1  0  0  0  0  1
 0  0  0  0  0  0  1  1  1  1  0  0  0
 0  0  1  0  1  0  1  0  0  0  0  1  0
 0  1  0  0  0  1  1  0  0  0  1  0  0
 0  0  0  1  1  1  0  0  0  1  0  0  0
 0  0  1  1  0  0  0  1  0  0  1  0  0
 0  1  0  1  0  0  0  0  1  0  0  1  0
 1  0  0  1  0  0  1  0  0  0  0  0  1
 1  0  0  0  0  1  0  1  0  0  0  1  0
 1  0  0  0  1  0  0  0  1  0  1  0  0
 1  1  1  0  0  0  0  0  0  1  0  0  0
 
julia> incidence_graph(3)
SimpleGraph{Int64} (n=26, m=52)

julia> H = incidence_hypergraph(3)
SimpleHypergraph{HVector{Mod{3}}} (n=13, m=13)
 ```

## Hilbert matrices

The file `hilbert.jl` contains code to generate Hilbert matrices.
```julia
julia> hilbert(4)
4×4 Array{Rational{Int64},2}:
 1//1  1//2  1//3  1//4
 1//2  1//3  1//4  1//5
 1//3  1//4  1//5  1//6
 1//4  1//5  1//6  1//7

julia> detx(ans)
1//6048000
 ```
