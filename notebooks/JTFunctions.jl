module JTFunctions

    ######
    using TaylorSeries, TaylorIntegration, Plots
    pyplot()
    import TaylorSeries: NumberNotSeries

    export evaluate_neighborhood, circle2, square2,
           ξmax, anal_vs_taylor2D, area_of_polygon,
           grid_ξmax, grid_FTLE, vectorField_plot,
           harmonic_oscillator!, simple_pendulum!, artificial_ode!
    ######

    ##### Functions & Indicators ####
    # Neighborhood Evaluation
    function evaluate_neighborhood{T<:NumberNotSeries}(δU::Function,ϕ::Taylor1{T},num_vals=100)
        τrange = linspace(0,1,num_vals)
        length(δU(τrange[1])) != 1 ? error("δU is not of the same dimension of the phase space.") : nothing

        return evaluate.(ϕ, δU.(τrange))
    end

    function evaluate_neighborhood{T<:NumberNotSeries}(δU::Function,ϕ::TaylorN{T},num_vals::Integer=100,ξ::Real=0.1)
        τrange = linspace(0,1,num_vals)
        length(δU(τrange[1])) != get_numvars() ? error("δU is not of the same dimension of the phase space.") : nothing
        !(δU(τrange[1]) ≈ δU(τrange[end])) ? error("δU is not a closed surface in the range given.") : nothing

        return evaluate.(ϕ, δU.(τrange,r=ξ))
        end

    function evaluate_neighborhood{T<:AbstractSeries}(δU::Function,ϕ::Array{T,1};num_vals::Integer=100,ξ::Real=0.1)
        return evaluate_neighborhood.(δU,ϕ,num_vals,ξ)
    end

    # Some Neighborhood Parametrization (this should deserve a special type Neighborhood)
    circle2(τ;r=0.1) = [r*cos(2π*τ),r*sin(2π*τ)]
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


    #### PLOTTING & GRIDS ####
    #No parallelization is done yet.

    #ξmax grid without the plot.. should the plot be included?
    #### This should be eliminated after PR is merged
    import TaylorSeries: get_variables
    get_variables(order) = [TaylorN(i,order=order) for i in 1:get_numvars()]
    ####
    function grid_ξmax{T<:Real}(eqs_diff!::Function,t0::T,tmax::T;
                        xlim::Tuple=(-1,1),ylim::Tuple=(-1,1),num_points::Integer=30,
                        order_jet::Integer=3,order_taylor::Integer=25,abstol=1e-10)
        @assert get_numvars() == 2 "ξmax must be 2-dimensional."
        δx,δy = get_variables(order_jet)

        xgrid= linspace(xlim...,num_points)
        ygrid= linspace(ylim...,num_points)
        heatgrid = zeros(length(xgrid),length(ygrid))

        #Parallelize here!
        for (i,x) in enumerate(xgrid)

            for (j,y) in enumerate(ygrid)
                q0TN = [x+δx,y+δy]
                _,ϕN = taylorinteg(eqs_diff!,q0TN,t0,tmax,order_taylor,abstol);
                heatgrid[i,j] = ξmax(ϕN)
            end

        end
        return xgrid, ygrid, heatgrid'
    end

    #FTLE scalar field using TaylorIntegration lyapunov expectrum.
    #TODO: FTLEs could also be made with MY method... (JT of order = 1 jets)
    function grid_FTLE{T<:Real}(eqs_diff!::Function,t0::T,tmax::T;
                            xlim::Tuple=(-1,1),ylim::Tuple=(-1,1),num_points::Integer=30,
                            order_taylor::Integer=25,abstol=1e-10)

        qgrid = linspace(xlim...,num_points)
        pgrid = linspace(ylim...,num_points)

        FTLE = zeros(num_points,num_points)

        #Parallelize here!
        for (i,q) in enumerate(qgrid)
            for (j,p) in enumerate(pgrid)

                x0 = [q,p]
                _,__,λ = liap_taylorinteg(eqs_diff!,x0,t0,tmax,order_taylor,abstol)
                FTLE[i,j] = maximum(λ[end,:])

            end
        end

        return qgrid, pgrid, FTLE'
    end

    #2D vector Field plotting for an TaylorIntegration defined eqs_motion!
    #See some options for beautifying (optional)
    function vectorField_plot(vecField!::Function,t0=0.0;xlim::Tuple=(-1,1),ylim::Tuple=(-1,1),
                              num_points::Int=15,arrowsize::Real=1.2,gridlog::Bool=false)
        X = linspace(xlim...,num_points)
        Y = linspace(ylim...,num_points)

        Δr = max(X.step.hi,Y.step.hi)
        sizeGrid = num_points^2

        #allocation
        vgrid = Vector{Tuple{Float64,Float64}}(sizeGrid)
        f_norm = Vector{Tuple{Float64,Float64}}(sizeGrid)
        normFactor = -Inf
        dx = zeros(Float64,2)

        for (i,x_) in enumerate(X)
            for (j,y_) in enumerate(Y)

                k = num_points*(i-1)+j
                vgrid[k] = (x_,y_)
                vecField!(t0,vgrid[k],dx)
                normFactor = max(normFactor, norm(dx))
                f_norm[k] = Tuple(arrowsize^2 * Δr * dx)
                if gridlog
                    f_norm[k] = sign.(f_norm[k]).*log.(norm.(f_norm[k]).+1)
                end

            end
        end

        if gridlog
            normFactor = log(normFactor+1)
        end

        f_norm = [f_norm[i]./normFactor for i in eachindex(f_norm)]

        #Plotting
        quiver(vgrid, quiver=f_norm)
        xlims!(xlim...)
        ylims!(ylim...)
    end


    #### ODE systems ####
    function harmonic_oscillator!(t,x,dx)
        dx[1] = -x[2]
        dx[2] = x[1]
        nothing
    end

    function simple_pendulum!(t,x,dx)
        θ,p = (x[1],x[2])
        dx[1] = p/(m*l^2)
        dx[2] = -m*g*l*sin(θ)
        nothing
    end
    m,g,l = 1.0,1.0,1.0

    function artificial_ode!(t,x,dx)
        dx[1] = 2.0*x[1]*x[2]
        dx[2] = -x[2]^2 + x[1]
        nothing
    end


end #module
