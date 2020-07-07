module LinearAlgebraX
using LinearAlgebra


# IntegerX is any sort of real or Gaussian integer
IntegerX = Union{S,Complex{S}} where S<:Integer

# RationalX is a Rational or Complex Rational based on integers 
RationalX = Union{Rational{S},Complex{Rational{S}}} where S<:Integer


function _recip(x::T) where T <: IntegerX
    return 1//x
end
_recip(x) = inv(x)


include("detx.jl")




end # module
