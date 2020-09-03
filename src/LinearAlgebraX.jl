module LinearAlgebraX
using LinearAlgebra, SimplePolynomials, Mods


# IntegerX is any sort of real or Gaussian integer
IntegerX = Union{S,Complex{S}} where S<:Integer

# RationalX is a Rational or Complex Rational based on integers
RationalX = Union{Rational{S},Complex{Rational{S}}} where S<:Integer

TypeX = Union{IntegerX, RationalX, AbstractMod}


function _recip(x::T) where T <: IntegerX
    return 1//x
end
_recip(x) = inv(x)

include("row_ops.jl")
include("detx.jl")
include("cofactor_det.jl")
include("eye.jl")
include("rrefx.jl")
include("invx.jl")
include("rankx.jl")
include("nullspacex.jl")
include("char_poly.jl")
include("homogeneous.jl")


end # module
