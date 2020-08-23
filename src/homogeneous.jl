import Base: show, eltype, length, size, getindex, ==, (*), iszero
import Base: Vector, Matrix, adjoint, hash

import LinearAlgebra: dot 
export HVector, HMatrix, dot, hash


"""
An `HVector` is a homogeneous vector; two `HVector`s are equal 
iff they are nonzero multiples of each other.

`HVector(a,b,c,...)` or `HVector([a,b,c,...])` creates a new `HVector `
with coordinates `(a,b,c,...)`.
"""
struct HVector{T}  
    data::Vector{T}
    function HVector(dat::Vector{T}) where T 
        data = dat .// one(T)
        S = typeof(sum(data))
        new{S}(_canonical(data))
    end 

end
HVector(x...) = HVector(collect(x))
HVector(x::HVector) = HVector(x.data)

"""
An `HMatrix` is a homogenous matrix; two of these are equal if they are
nonzero multiples of each other. 

Given a two-dimensional array `A`, `HMatrix(A)` creates a new homogeneous 
matrix. 
"""
struct HMatrix{T}
    data::Matrix{T}
    function HMatrix(dat::AbstractMatrix{T}) where T 
        data = dat .// one(T)
        S = typeof(sum(data))
        M = new{S}(_canonical(data))
        return M 
    end 
end
HMatrix(x::HMatrix) = HMatrix(x.data)


function _canonical(x)
    if length(x) == 0 || iszero(x)
        return x 
    end 
    idx = findlast(x .!= 0)
    a = x[idx]
    return x .// a 
end 

HObject = Union{HVector,HMatrix}


eltype(x::HObject) = eltype(x.data)
length(x::HObject) = length(x.data)
size(x::HObject) = size(x.data)
getindex(x::HObject,k...) = x.data[k...]
(==)(x::HObject,y::HObject) = x.data == y.data
iszero(x::HObject) = iszero(x.data)


"""
The `dot` produce of two `HVector`s is either `0` 
(if they are orthogonal) or `1` (otherwise).
"""
function dot(x::HVector, y::HVector)::Int
    d = dot(x.data,y.data)
    if d==0
        return 0 
    end
    return 1
end 


function show(io::IO, x::HVector)
    print(io,"HVector(")
    n = length(x.data)
    for j=1:n-1
        print(io,"$(x.data[j]), ")
    end 
    print(io,x.data[end],")")
end

function show(io::IO, x::HMatrix)
    print(io,"HMatrix: ")
    print(io,x.data)
end

function (*)(A::AbstractArray{T,2}, v::HVector) where T 
    Av = A*v.data 
    return HVector(Av)
end 

function (*)(A::AbstractArray{T,2}, B::HMatrix) where T 
    AB = A*B.data 
    return HMatrix(AB)
end 

(*)(A::HMatrix, v::HVector) = HVector(A.data * v.data)
(*)(A::HMatrix, B::HMatrix) = HMatrix(A.data * B.data)

adjoint(A::HMatrix) = HMatrix(adjoint(A.data))

hash(x::HObject,h::UInt=UInt(0)) = hash(x.data,h)

Vector(x::HVector) = copy(x.data)
Matrix(x::HMatrix) = copy(x.data)