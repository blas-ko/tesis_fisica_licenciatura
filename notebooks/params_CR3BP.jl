### CR3BP Functions ###
V_pseudo(r) = _V_pseudo(r,μ)
CR3BP_ode!(t, x, dx) = _CR3BP_ode!(t, x, dx, μ)
C_Jacobi(r,v) = _C_Jacobi(r,v,μ)
vy(r, vx; E=C_L1) = _vy(r, vx, μ, E=E)
### CR3BP params ###
L_1  =  1 - (μ/3)^(1/3)
L_2  =  1 + (μ/3)^(1/3)
L_3  =  -(1 + (5μ/12))
L_x4 = 1/2 * (1 - 2μ)
L_x5 = 1/2 * (1 - 2μ)

C_L1 = V_pseudo([L_1,0.0])
C_L2 = V_pseudo([L_2,0.0])
C_L3 = V_pseudo([L_3,0.0])
C_L4 = V_pseudo([L_x4,L_y4])
C_L5 = V_pseudo([L_x5,L_y5])

println("L_x4: $(round(L_x4,3))")
println("C_L1: $(round(C_L1,3))")


### Integration Parameters ###
abstol = 1e-7
ord_taylor = 20
ord_jet = 4


println("\nTaylorIntegration parameters:
abstol = $abstol
ord_taylor = $ord_taylor
ord_jet = $ord_jet")

### Independent configuration variables ###
δq1, δq2 = set_variables("x y",order=ord_jet)

println("\nIndependent Variables:
δq1: $δq1
δq2: $δq2")

println("\nμ: $μ")
