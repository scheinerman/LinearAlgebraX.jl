# Extras Folder

## Projective planes

The file `projective.jl` contains code for creating structures related 
to finite projective planes of prime order.

* `generate_points(p)` generates a list of the homogeneous coordinates of the 
points (lines) of a finite projective plane of prime order `p`.
* `incidence_matrix(p)` generates the point-line incidence matrix of a 
projective plane of order `p`.
* `incidence_graph(p)` generates the bipartite point-line incidence graph 
of a finite projective plane of order `p`.
```julia
julia> include("extras/projective.jl")
incidence_graph

julia> A = incidence_matrix(3)
13×13 Array{Int64,2}:
 1  1  1  1  1  1  1  1  1  0  0  0  0
 1  1  0  1  1  0  1  1  0  1  1  1  0
 1  0  1  1  0  1  1  0  1  1  1  1  0
 1  1  1  1  1  1  0  0  0  0  1  1  1
 1  1  0  1  0  1  0  1  1  1  1  0  1
 1  0  1  1  1  0  0  1  1  1  0  1  1
 1  1  1  0  0  0  1  1  1  0  1  1  1
 1  1  0  0  1  1  1  0  1  1  0  1  1
 1  0  1  0  1  1  1  1  0  1  1  0  1
 0  1  1  0  1  1  0  1  1  1  1  1  0
 0  1  1  1  1  0  1  0  1  1  1  0  1
 0  1  1  1  0  1  1  1  0  1  0  1  1
 0  0  0  1  1  1  1  1  1  0  1  1  1

julia> G = incidence_graph(3)
SimpleGraph{Int64} (n=26, m=117)
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
 ```