module JTFunctions

    using TaylorSeries, TaylorIntegration
    import TaylorSeries: NumberNotSeries

    export evaluate_neighborhood, circle2, square2,
           ξmax, anal_vs_taylor2D, area_of_polygon

    # Neighborhood Evaluation
    function evaluate_neighborhood{T<:NumberNotSeries}(δU::Function,ϕ::Taylor1{T},num_vals=100)
        τrange = linspace(0,1,num_vals)
        length(δU(τrange[1])) != 1 ? error("δU is not of the same dimension of the phase space.") : nothing

        return evaluate.(ϕ, δU.(τrange))
    end

    function evaluate_neighborhood{T<:NumberNotSeries}(δU::Function,ϕ::TaylorN{T},num_vals=100)
        τrange = linspace(0,1,num_vals)
        length(δU(τrange[1])) != get_numvars() ? error("δU is not of the same dimension of the phase space.") : nothing
        !(δU(τrange[1]) ≈ δU(τrange[end])) ? error("δU is not a closed surface in the range given.") : nothing

        return evaluate.(ϕ, δU.(τrange))
        end

    function evaluate_neighborhood{T<:AbstractSeries}(δU::Function,ϕ::Array{T,1},num_vals=100)
        return evaluate_neighborhood.(δU,ϕ,num_vals)
    end

    # Some Neighborhood Parametrization (this should deserve a special type Neighborhood)
    circle2(τ,r=0.1) = [r*cos(2π*τ),r*sin(2π*τ)]
    square2(τ) = nothing

    # Area of Neighborhood
    function area_of_polygon(pjet,qjet)
        #points_per_polygon = length(pjet[1])
        num_time_steps = length(pjet)

        area_per_polygon = zeros(eltype(pjet[1]), num_time_steps)
        for i in 1:num_time_steps

            p_aux = copy(pjet[i])

            p_end = shift!(p_aux)
            push!(p_aux,p_end)

            A_1 = sum(p_aux.*qjet[i])

            q_aux = copy(qjet[i])

            q_end = shift!(q_aux)
            push!(q_aux,q_end)

            A_2 = sum(q_aux.*pjet[i])

            area_per_polygon[i] = (A_1 - A_2)/2.0
        end
        return area_per_polygon
    end

    # Maximum Box Size
    function ξmax(ϕ;ϵ=10e-6)
        dof = size(ϕ)[2]
        ord = ϕ[end,1].order
        ξ = Inf

        for m in 1:dof
            a_k = ϕ[end,m].coeffs[ord+1].coeffs

            #In case the last term is zero.
            k = 1
            while a_k == zeros(a_k)
                a_k = ϕ[end,m].coeffs[ord+1-k].coeffs
                k += 1
            end

            for a_jk in a_k
                ξ = min(ξ,abs(ϵ/a_jk)^(1./ord) )
            end
        end

        return ξ
    end


    # Compare taylorinteg solution `x,y` vs analytical ones `xa`
    anal_vs_taylor2D(x,y,xa) = sqrt.((x - xa[:,1]).^2 + (y - xa[:,2]).^2)

end #module
