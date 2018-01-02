using TaylorSeries
import TaylorSeries: evaluate 

#Taylor1 matrix evaluation 
function evaluate(A::Array{Taylor1{T},2}, δt::S) where {T<:Number, S<:Number}
    R = promote_type(T,S)
    return evaluate(convert(Array{Taylor1{R},2},A), convert(R,δt))
end
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
evaluate(A::Array{Taylor1{T},2}) where {T<:Number} = evaluate(A,zero(T))

(p::Array{Taylor1{T},2})(x) where {T<:Number} = evaluate(p, x)

(p::Array{Taylor1{T},2})() where {T<:Number} = evaluate.(p)

#promotion of TaylorN Vector evaluation
function evaluate(A::Array{TaylorN{T},1}, δx::Vector{S}) where {T<:Number, S<:Number}
    R = promote_type(T,S)
    return evaluate(convert(Array{TaylorN{R},1},A), convert(Vector{R},δx))
end

#TaylorN matrix evaluation
function evaluate(A::Array{TaylorN{T},2}, δx::Vector{S}) where {T<:Number, S<:Number}
    R = promote_type(T,S)
    return evaluate(convert(Array{TaylorN{R},2},A), convert(Vector{R},δx))
end
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
evaluate(A::Array{TaylorN{T},2}) where {T<:Number} = evaluate.(A)

(p::Array{TaylorN{T},2})(x) where {T<:Number} = evaluate(p, x)

(p::Array{TaylorN{T},2})() where {T<:Number} = evaluate.(p)