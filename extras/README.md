# Extras Folder

 The `Extras` folder contains `projective.jl` for creating structures related 
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
