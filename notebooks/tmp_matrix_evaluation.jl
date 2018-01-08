using TaylorSeries
import TaylorSeries: evaluate

function evaluate(x::Union{Array{Taylor1{T},1}, SubArray{Taylor1{T},1}}, δt::S) where {T<:Number, S<:Number}
    R = promote_type(T,S)
    return evaluate(convert(Array{Taylor1{R},1},x), convert(R,δt))
end
# function evaluate(x::SubArray{Taylor1{T},1}, δt::S) where {T<:Number, S<:Number}
#     R = promote_type(T,S)
#     return evaluate(convert(SubArray{Taylor1{R},1},x), convert(R,δt))
# end
# function evaluate(x::SubArray{Taylor1{T},1}, δt::T) where {T<:Number}
#     xnew = Array{T}( length(x) )
#     evaluate!(x, δt, xnew)
#     return xnew
# end

#Taylor1 matrix evaluation
function evaluate(A::Union{Array{Taylor1{T},2}, SubArray{Taylor1{T},2}}, δt::S) where {T<:Number, S<:Number}
    R = promote_type(T,S)
    return evaluate(convert(Array{Taylor1{R},2},A), convert(R,δt))
end
# function evaluate(A::SubArray{Taylor1{T},2}, δt::S) where {T<:Number, S<:Number}
#     R = promote_type(T,S)
#     return evaluate(convert(SubArray{Taylor1{R},2},A), convert(R,δt))
# end

function evaluate(A::Array{Taylor1{T},2}, δt::T) where {T<:Number}
    n,m = size(A)
    Anew = Array{T}( n,m )
    xnew = Array{T}( n )

    for i in 1:m
        evaluate!(A[:,i], δt, xnew)
        Anew[:,i] = xnew
    end

    return Anew
end
# function evaluate(A::SubArray{Taylor1{T},2}, δt::T) where {T<:Number}
#     n,m = size(A)
#     Anew = Array{T}( n,m )
#     xnew = Array{T}( n )
#
#     for i in 1:m
#         evaluate!(A[:,i], δt, xnew)
#         Anew[:,i] = xnew
#     end
#
#     return Anew
# end
evaluate(A::Array{Taylor1{T},2}) where {T<:Number} = evaluate.(A)
evaluate(A::SubArray{Taylor1{T},2}) where {T<:Number} = evaluate.(A)

(p::Array{Taylor1{T},2})(x) where {T<:Number} = evaluate(p, x)
(p::SubArray{Taylor1{T},2})(x) where {T<:Number} = evaluate(p, x)
(p::Array{Taylor1{T},2})() where {T<:Number} = evaluate.(p)
(p::SubArray{Taylor1{T},2})() where {T<:Number} = evaluate.(p)

#promotion of TaylorN Vector evaluation
function evaluate(A::Union{Array{TaylorN{T},1}, SubArray{TaylorN{T},1}}, δx::Vector{S}) where {T<:Number, S<:Number}
    R = promote_type(T,S)
    return evaluate(convert(Array{TaylorN{R},1},A), convert(Vector{R},δx))
end
# function evaluate(A::SubArray{TaylorN{T},1}, δx::Vector{S}) where {T<:Number, S<:Number}
#     R = promote_type(T,S)
#     return evaluate(convert(SubArray{TaylorN{R},1},A), convert(Vector{R},δx))
# end

# function evaluate(x::SubArray{TaylorN{T},1}, δx::Array{T,1}) where {T<:Number}
#     x0 = Array{T}( length(x) )
#     evaluate!( x, δx, x0 )
#     return x0
# end


#TaylorN matrix evaluation
function evaluate(A::Union{Array{TaylorN{T},2},SubArray{TaylorN{T},2}} , δx::Vector{S}) where {T<:Number, S<:Number}
    R = promote_type(T,S)
    return evaluate(convert(Array{TaylorN{R},2},A), convert(Vector{R},δx))
end
# function evaluate(A::SubArray{TaylorN{T},2}, δx::Vector{S}) where {T<:Number, S<:Number}
#     R = promote_type(T,S)
#     return evaluate(convert(SubArray{TaylorN{R},2},A), convert(Vector{R},δx))
# end

function evaluate(A::Array{TaylorN{T},2}, δx::Vector{T}) where {T<:Number}
    n,m = size(A)
    Anew = Array{T}( n,m )
    xnew = Array{T}( n )

    for i in 1:m
        evaluate!(A[:,i], δx, xnew)
        Anew[:,i] = xnew
    end

    return Anew
end
# function evaluate(A::SubArray{TaylorN{T},2}, δx::Vector{T}) where {T<:Number}
#     n,m = size(A)
#     Anew = Array{T}( n,m )
#     xnew = Array{T}( n )
#
#     for i in 1:m
#         evaluate!(A[:,i], δx, xnew)
#         Anew[:,i] = xnew
#     end
#
#     return Anew
# end

evaluate(A::Array{TaylorN{T},2}) where {T<:Number} = evaluate.(A)
evaluate(A::SubArray{TaylorN{T},2}) where {T<:Number} = evaluate.(A)

(p::Array{TaylorN{T},2})(x) where {T<:Number} = evaluate(p, x)
(p::SubArray{TaylorN{T},2})(x) where {T<:Number} = evaluate(p, x)
(p::Array{TaylorN{T},2})() where {T<:Number} = evaluate.(p)
(p::SubArray{TaylorN{T},2})() where {T<:Number} = evaluate.(p)
