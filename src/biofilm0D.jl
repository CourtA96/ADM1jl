export petersenmatrixtranspose_biofilmAdditions
function petersenmatrixtranspose_biofilmAdditions()
    
   # Attachemnt and detachment
   P2_bio = [ # columns of the Petersen matrix
   13,13, # suspended
   14,14,
   15,15,
   16,16,
   17,17,
   18,18,
   19,19,
   20,20,
   21,21,
   22,22,
   23,23,
   24,24,
   36,36, # biofilm
   37,37,
   38,38,
   39,39,
   40,40,
   41,41,
   42,42,
   43,43,
   44,44,
   45,45,
   46,46,
   47,47,
   ]

   P1_bio = [ # rows of the Petersen matrix
   1,13, #1 attach, detach
   2,14, #2
   3,15, #3
   4,16, #4
   5,17, #5
   6,18, #6
   7,19, #7
   8,20, #8
   9,21, #9
   10,22, #10
   11,23, #11
   12,24, #12
   1,13, #36 attach, detach
   2,14, #37
   3,15, #38
   4,16, #39
   5,17, #40
   6,18, #41
   7,19, #42
   8,20, #43
   9,21, #44
   10,22, #45
   11,23, #46
   12,24, #47
   ]

   Q_bio = [ # attachment, detachment
   -1,1, #1 suspended
   -1,1, #2
   -1,1, #3
   -1,1, #4
   -1,1, #5
   -1,1, #6
   -1,1, #7
   -1,1, #8
   -1,1, #9
   -1,1, #10
   -1,1, #11
   -1,1, #12
   1,-1, #36 biofilm
   1,-1, #37
   1,-1, #38
   1,-1, #39
   1,-1, #40
   1,-1, #41
   1,-1, #42
   1,-1, #43
   1,-1, #44
   1,-1, #45
   1,-1, #46
   1,-1, #47
   ]

   PM_bio = sparse(P2_bio, P1_bio, Q_bio) # return transpose

   return Matrix(PM_bio)
end

export biofilmRates
function biofilmRates(u,AD)
    total_suspend = sum(u[13:24])
    total_biofilm = sum(u[36:47])
    ka = AD[1] # attachment coefficient 0.3
    kd = AD[2] # detachment coefficient 0.01

    rrates = Array{Real}(undef,24)
    rrates[1:12] = ka*total_suspend*u[13:24] # attachment rates
    rrates[13:24] = kd*total_biofilm*u[36:47] # detachment rates

    return rrates
    
end

export biofilmRHSfun
"""
    biofilmRHSfun(du,u,p,t)

Return the right-hand side of the system of 0D biofilm ODEs, this is an in-place function.

# Arguments
- `du::Vector`: the rate change of the state vector (required since the function is defined in-place).
- `u::Vector`: the state vector.
- `p::Vector`: all of the model parameters.
- `t`: the timestep, usually a Float64.
"""
function biofilmRHSfun(du,u,p,t)
   # Parameters:

   RP = p[1:6] # reactor parameters
   BP = p[7:51] # biochemical paramters
   SP = p[52:69] # stoichiometric parameters
   CC = p[70:87] # carbon balance parameters
   PhP = p[88:107] # physiochemical parameters
   IV_temp = p[108:142] # inflow vector

   IV = [IV_temp; zeros(12)]

   t1 = p[143] # start time
   tMax = p[144] # maximum time

   AD = p[145:146] # attach and detach params

   # compute pressures
   pressures = pressureOfGasses(u,PhP,RP)

   # compute transport matrix
   TM = transportmatrix_definition(RP,pressures)

   # create Petersen Matrix
   PM_temp = petersenmatrixtranspose_definition(RP,BP,SP,CC)
   PM_bio = petersenmatrixtranspose_biofilmAdditions()

   PM = [PM_temp PM_bio]

   # compute reaction rates
   rr_temp = reactionrates(BP,RP,PhP,pressures,u,48)
   rr_bio = biofilmRates(u,AD)

   # rr = [rr_temp;rr_bio]
   rr = [rr_temp; rr_bio]

   # compute inflow/outflow rate
   liquidFlow = RP[6]/RP[4]

   t2 = time()

   if t2-t1 > tMax
      errStr = string("Took longer than ",tMax, " seconds.")
      error(errStr)
   end

   du .= TM*u + PM*rr + IV*liquidFlow

end

export ExampleBiofilm
"""
    function ExampleBiofilm(tspan::Tuple,u0::Vector,IV::Vector; <keyword arguments>)

Compute the solution for the given timespan, `tspan`; initial condition, `u0`; and inflow vector `IV`.

Also return the time (in seconds) the solution took to compute.

# Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.

# Examples
```julia-repl
julia> u0 = ADM1code.InitialConditions();

julia> IV = ADM1code.inflowvector_definition();

julia> sol, tSol = ADM1code.ExampleSol((0.0,200.0),u0,IV);

julia> sol
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 146-element Vector{Float64}:
[...]

u: 146-element Vector{Vector{Float64}}:
[...]

julia> tSol
0.3854937
```
"""
function ExampleBiofilm(tspan::Tuple,u0::Vector,IV::Vector; tols=1e-4,tMax = 300.0)

   # Combine all of the parameter vectors into one vector (parm) for the solver.

   RP = reactorParameterDefinition() # length = 6
   BP = biochemicalparameter_definition() # length = 45
   SP = stoichiometricparameter_definition() # length = 18
   CC = carbonContent_definition() # length = 18
   PhP = physiochemicalParameterDefinition(RP) # length = 20
   #IV = inflowvector_definition() # length = 35

   parm = zeros(146)
   parm[1:6] = RP
   parm[7:51] = BP
   parm[52:69] = SP
   parm[70:87] = CC
   parm[88:107] = PhP
   parm[108:142] = IV
   parm[143] = time() # start time
   parm[144] = tMax # maximum time
   parm[145:146] = [0.3,0.01]

   prob=ODEProblem(biofilmRHSfun,u0,tspan,parm)

   global sol = "not defined"

   # try
      global t = @timed sol = solve(prob,alg=Rosenbrock23())
   # catch e
   #    printstyled(stderr,"\nERROR: ", bold=true, color=:red)
   #    printstyled(stderr,sprint(showerror,e), color=:light_red)
   #    println(stderr)
   # end

   # isSolved = (sol != "not defined")

   # if isSolved == true
      return sol,t[2]
   # else
   #    # return nothing # normally the method should return nothing since the problem was not solved
   #    return ["tMax reached" for i in 1:47],"tMax reached" # only use this for testing
   # end

end

export Biofilm
"""
    function Biofilm(tspan::Tuple,u0::Vector,IV::Vector; <keyword arguments>)

Compute the solution for the given timespan, `tspan`; initial condition, `u0`; and inflow vector `IV`.

Also return the time (in seconds) the solution took to compute.

# Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.

# Examples
```julia-repl
julia> u0 = ADM1code.InitialConditions();

julia> IV = ADM1code.inflowvector_definition();

julia> sol, tSol = ADM1code.ExampleSol((0.0,200.0),u0,IV);

julia> sol
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 146-element Vector{Float64}:
[...]

u: 146-element Vector{Vector{Float64}}:
[...]

julia> tSol
0.3854937
```
"""
function Biofilm(tspan::Tuple,u0::Vector,IV::Vector; tols=1e-4,tMax = 300.0)

   # Combine all of the parameter vectors into one vector (parm) for the solver.

   # Combine all of the parameter vectors into one vector (parm) for the solver.

   p = CSV.File("model_parameters.csv")
   # values indexed by column header, ie: params_temp.R  returns 0.083145
   RP = [
   p.T_ad, p.T_base,
   p.P_atm,
   p.V_liq, p.V_gas,
   p.Q_ad
   ]

   BP = [p.k_dis, p.k_hyd_ch, p.k_hyd_pr, p.k_hyd_li, p.tresX,
         p.k_dec_all, p.K_S_IN,
         p.pH_UL_aa, p.pH_LL_aa,
         p.k_m_su, p.K_S_su, p.Y_su,
         p.k_m_aa, p.K_S_aa, p.Y_aa,
         p.k_m_fa, p.K_S_fa, p.Y_fa, p.K_I_h2_fa,
         p.k_m_c4, p.K_S_c4, p.Y_c4, p.K_I_h2_c4,
         p.k_m_pro, p.K_S_pro, p.Y_pro, p.K_I_h2_pro,
         p.k_m_ac, p.K_S_ac, p.Y_ac, p.pH_UL_ac, p.pH_LL_ac, p.K_I_nh3,
         p.k_m_h2, p.K_S_h2,
         p.Y_h2, p.pH_UL_h2, p.pH_LL_h2,
         p.k_dec_X_su, p.k_dec_X_aa, p.k_dec_X_fa, p.k_dec_X_c4, p.k_dec_X_pro, p.k_dec_X_ac, p.k_dec_X_h2
         ]

   SP =  [
         p.f_sI_xc, p.f_xI_xc, p.f_ch_xc, p.f_pr_xc, p.f_li_xc,
         p.N_xc, p.N_I, p.f_fa_li,
         p.f_h2_su, p.f_bu_su, p.f_pro_su, p.f_ac_su, p.f_h2_aa,
         p.N_aa, p.f_va_aa, p.f_bu_aa, p.f_pro_aa, p.f_ac_aa
         ] # Note that N_xs and N_I values differ between python and julia


   C_IC = 1.0
   C_IN = 0.0
   C_h2 = 0.0
   CC = [
         p.C_su, p.C_aa, p.C_fa, p.C_va, p.C_bu, p.C_pro, p.C_ac, C_h2, p.C_ch4, C_IC, C_IN, p.C_sI,
         p.C_xc, p.C_ch, p.C_pr, p.C_li, p.C_bac, p.C_xI
         ]

   PhP = [
         p.R,
         p.K_h2o, p.k_AB_va, p.k_AB_bu, p.k_AB_pro, p.k_AB_ac, p.k_AB_co2, p.k_AB_IN,
         p.K_a_va, p.K_a_bu, p.K_a_pro, p.K_a_ac, p.K_a_co2, p.K_a_IN,
         p.p_gas_h2o, p.k_p,
         p.k_L_a, p.K_H_co2, p.K_H_ch4, p.K_H_h2
         ]

   AD = [
      p.ka, p.kd
      ]

   RP = [i[1] for i in RP]
   BP = [i[1] for i in BP]
   SP = [i[1] for i in SP]
   CC = [i[1] for i in CC]
   PhP = [i[1] for i in PhP]
   AD = [i[1] for i in AD]

   parm = zeros(146)
   parm[1:6] = RP
   parm[7:51] = BP
   parm[52:69] = SP
   parm[70:87] = CC
   parm[88:107] = PhP
   parm[108:142] = IV
   parm[143] = time() # start time
   parm[144] = tMax # maximum time
   parm[145:146] = AD

   prob=ODEProblem(biofilmRHSfun,u0,tspan,parm)

   global sol = "not defined"

   # try
      global t = @timed sol = solve(prob,alg=Rosenbrock23())
   # catch e
   #    printstyled(stderr,"\nERROR: ", bold=true, color=:red)
   #    printstyled(stderr,sprint(showerror,e), color=:light_red)
   #    println(stderr)
   # end

   # isSolved = (sol != "not defined")

   # if isSolved == true
      return sol,t[2]
   # else
   #    # return nothing # normally the method should return nothing since the problem was not solved
   #    return ["tMax reached" for i in 1:47],"tMax reached" # only use this for testing
   # end

end