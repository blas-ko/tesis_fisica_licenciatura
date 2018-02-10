module ThreeBody

    ######
    using TaylorSeries, TaylorIntegration

    export M_1, M_2, μ, mu,
        L_1,L_2,L_3,L_x4,L_y4,L_x5,L_y5,
        L1,L2,L3,Lx4,Lx5,
        C_L1,C_L2,C_L3,C_L4,C_L5,
        mass_moon, mass_sun, mass_earth,
        _V_pseudo, _CR3BP_ode!, _C_Jacobi, _vy,
        collision_approach, JT_accuracy

    ######

    ### Earth, Moon, Sun masses ###
    const mass_moon  = 7.347673e22 #kg
    const mass_earth = 5.9722e24 #kg
    const mass_sun = 1.9855e30 #kg
    # mass_satelite ≈ 1e3 kg
    ### Earth, Moon, Sun distances ###

    ### Mass parameters (with distance normalized to 1)
    M_1 = mass_earth
    M_2 = mass_moon
    mu(M_1,M_2) = M_2/(M_1 + M_2)
    μ = mu(M_1,M_2)



    ### Lagrangian Points ###
    # L1, L2, L3 can be obtained numerically more precisely
    """
    \$ 1 - ( μ /3 )^{1/3} \$
    """
    L1(μ)  =  1 - (μ/3)^(1/3)
    """
    \$ 1 + ( μ /3 )^{1/3} \$
    """
    L2(μ)  =  1 + (μ/3)^(1/3)
    """
    \$ - 1 - 5μ /12  \$
    """
    L3(μ)  =  -(1 + (5μ/12))
    """
    \$  (1 - 2μ)/2  \$
    """
    Lx4(μ) = 1/2 * (1 - 2μ)
    """
    \$  (1 - 2μ)/2  \$
    """
    Lx5(μ) = 1/2 * (1 - 2μ)

    L_1  = L1(μ)
    L_2  = L2(μ)
    L_3  = L3(μ)
    L_x4 = Lx4(μ)
    L_x5 = Lx5(μ)
    L_y4 = √3/2
    L_y5 = -√3/2

    ####  Potentials ###
    V_G(r, μ) = -(1-μ)/norm(r-[-μ,0.0]) - μ/norm(r-[1-μ,0.0])
    V_Ω(r) = -norm(r)^2/2

    """
    Returns the no coriolis potential for the CR3BP at value `r`.
    """
    _V_pseudo(r, μ) = V_G(r,μ) + V_Ω(r)


    # V_pseudo(r) = _V_pseudo(r,μ)

    ### Jacobi Constants ###
    C_L1 = _V_pseudo([L_1,0.0],μ)
    C_L2 = _V_pseudo([L_2,0.0],μ)
    C_L3 = _V_pseudo([L_3,0.0],μ)
    C_L4 = _V_pseudo([L_x4,L_y4],μ)
    C_L5 = _V_pseudo([L_x5,L_y5],μ);

    ### Eqs of motion ###
    """
    Equations of motion for CR3BP ready for integrating via `TaylorIntegration`
    """
    function _CR3BP_ode!(t, x, dx, μ)

        β = 1.0-μ
        r_13_cubed = ( (x[1]+μ)^2 + x[2]^2 )^(3/2)
        r_23_cubed = ( (x[1]-β)^2 + x[2]^2 )^(3/2)

        dx[1] = x[3]
        dx[2] = x[4]
        dx[3] = -( β*(x[1]+μ)/r_13_cubed + μ*(x[1]-β)/r_23_cubed ) + x[1] + 2x[4]
        dx[4] = -(β/r_13_cubed + μ/r_23_cubed - 1)*x[2] - 2x[3]

        nothing
    end


    # CR3BP_ode!(t, x, dx) = _CR3BP_ode!(t, x, dx, μ)


    """
    Distance between two 2-dimn solutions.
    """
    function JT_accuracy(ϕ,ϕ_JT)
        dims = size(ϕ)[2]
        diffs = zero(ϕ[:,1])
        for dim in 1:dims
            diffs += (ϕ[:,dim] - ϕ_JT[:,dim]).^2
        end
        return sqrt.(diffs)
    end

    """
    Returns the indexes for which two trajectories differ in less than `d_min`.
    """
    function collision_approach(ϕ1, ϕ2; d_min=0)
        d_ϕ1ϕ2 = JT_accuracy(ϕ1[:,1:2], ϕ2[:,1:2])
        if d_min == 0
            d_min = minimum(d_ϕ1ϕ2)
        end
        #approached_distance = Vector{eltype(d_A1_A2)}()
        #approached_time = Vector{eltype(tv)}()
        approached_index = Vector{Int64}()
        for i in eachindex(d_ϕ1ϕ2)
            if d_ϕ1ϕ2[i] <= d_min
                #push!(approached_distance, d_ϕ1ϕ2[i])
                #push!(approached_time, tv[i])
                push!(approached_index, i)
            end
        end

        return approached_index
    end

    ### Jacobi Constant ###
    """
    Computes Jacobi Energy \$ V(r) + v^2/2 \$.
    """
    _C_Jacobi(r,v,μ) = _V_pseudo(r,μ) + (v⋅v)/2.0

    ### implicit y-velocity ###
    """
    Obtains y-velocity for initial conditions to be at energy level `E`
    """
    _vy(r, vx, μ; E=C_L1) = √( 2*( E - _V_pseudo(r,μ) ) - vx^2 )

    println("M_1: $M_1 , M_2: $M_2")
    println("L_1: $(round(L_1,3)), L_2: $(round(L_2,3)), L_3: $(round(L_1,3)), L_4x: $(round(L_x4,3)), L_4y: $(round(L_y4,3))")
end
