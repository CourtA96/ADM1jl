export RHSfunInflowVaried
"""
    RHSfunInflowVaried(du,u,p,t)

Return the right-hand side of the system of ODEs, this is an in-place function.

# Arguments
- `du::Vector`: the rate change of the state vector (required since the function is defined in-place).
- `u::Vector`: the state vector.
- `p::Vector`: all of the model parameters.
- `t`: the timestep, usually a Float64.
"""
function RHSfunInflowVaried(du::Vector,u::Vector,p::Vector,t)
   # Parameters:

   RP = p[1:6] # reactor parameters
   BP = p[7:51] # biochemical paramters
   SP = p[52:69] # stoichiometric parameters
   CC = p[70:87] # carbon balance parameters
   PhP = p[88:107] # physiochemical parameters
   IV = [inflowFunctions[i][t] for i in 1:35] # inflow is the value of the interpolation of the previous cstr at time t

   t1 = p[108] # start time
   tMax = p[109] # maximum time

   # compute pressures
   pressures = pressureOfGasses(u,PhP,RP)

   # compute transport matrix
   TM = transportmatrix_definition(RP,pressures)

   # create Petersen Matrix
   PM = petersenmatrixtranspose_definition(RP,BP,SP,CC)

   # compute reaction rates
   rr = reactionrates(BP,RP,PhP,pressures,u,size(PM)[2])

   t2 = time()

   if t2-t1 > tMax
      errStr = string("Took longer than ",tMax, " seconds.")
      error(errStr)
   end

   du .= TM*u + PM*rr - TM*IV

end

"""
##########################################################
# MAIN PROGRAM ###########################################
##########################################################
"""

export MultiChamberSolutionExample
"""
    function MultiChamberSolutionExample(tspan::Tuple,u0::Tuple,IV::Vector,nChambers::Int64; <keyword arguments>)

Compute the solution for a system of `nChambers` connected CSTRs with `u0` initial conditions. The 
inflow of the first CSTR is given by `IV`, the outflow of the first CSTR becomes the inflow of second CSTR 
and so on.

# Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.

# Examples
```julia-repl
julia> u0 = ADM1code.initialConditions();

julia> IV = ADM1code.inflowvector_definition();

julia> sol, tSol = ADM1code.ExampleMultiChamberSol((0.0,200.0),u0,IV,3);

julia> sol
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 146-element Vector{Float64}:
[...]

u: 146-element Vector{Vector{Float64}}:
[...]

julia> tSol
0.3854937
"""
function MultiChamberSolutionExample(tspan::Tuple,u0::Tuple,IV::Vector,nChambers::Int64;alg = Rodas4P(), tols=1e-4,tMax = 300.0)
    
    sols = Vector{Any}(undef,nChambers)
    tSols = Vector{Any}(undef,nChambers)

    sols[1],tSols[1] = ExampleSol(tspan,u0[1],IV,alg=alg,tols=tols,tMax=tMax)

    println("Finished Chamber 1")

    for i in 2:nChambers
        if typeof(sols[i-1]) != Vector{String}
            sols[i],tSols[i] = ExampleMultiChamberSol(tspan,u0[i],sols[i-1],alg=alg,tols=tols,tMax=tMax)
            println("Finished Chamber ",i)
        else
            print("Could not complete chamber because previous chamber errored out.")
            sols[i] = ["tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached"]
        end
    end

    return sols,sum(tSols)

end

export MultiChamberSolution
"""
    function MultiChamberSolution(tspan::Tuple,u0::Tuple,IV::Vector,nChambers::Int64; <keyword arguments>)

Compute the solution for a system of `nChambers` connected CSTRs with `u0` initial conditions. The 
inflow of the first CSTR is given by `IV`, the outflow of the first CSTR becomes the inflow of second CSTR 
and so on.

This function reads in the parameter values from a .csv file.
The names of the .csv files that contain the parameter values should be "model_parameters.csv", "model_parameters2.csv", 
"model_parameters3.csv", ... and so on.

# Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.
- `saveAt = []`: specific times to save the solution. If given a number `n`, the solver will save the solution every `n` timesteps

# Examples
```julia-repl
julia> u0 = initialConditions();

julia> IV = inflowvector_definition();

julia> sols = MultiChamberSolution((0.0,200.0),(u0,u0,u0),IV,3);
Finished Chamber 1
Finished Chamber 2
Finished Chamber 3

julia> sols[1]
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 115-element Vector{Float64}:
[...]

u: 115-element Vector{Vector{Float64}}:
[...]

julia> sols[2]
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 115-element Vector{Float64}:
[...]

u: 115-element Vector{Vector{Float64}}:
[...]

julia> sols[3]
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 115-element Vector{Float64}:
[...]

u: 115-element Vector{Vector{Float64}}:
[...]
```
"""
function MultiChamberSolution(tspan::Tuple,u0::Tuple,IV::Vector,nChambers::Int64;alg = Rodas4P(), tols=1e-4,tMax = 300.0,saveAt=[])
    
    sols = Vector{Any}(undef,nChambers)

    sols[1] = ADM1sol(tspan,u0[1],IV,alg=alg,tols=tols,tMax=tMax,saveAt=saveAt)[1]

    println("Finished Chamber 1")

    for i in 2:nChambers
        if typeof(sols[i-1]) != Vector{String}
            filename = string("model_parameters",i,".csv")
            sols[i] = ADM1MultiChamberSol(tspan,u0[i],sols[i-1],paramFilename = filename,alg=alg,tols=tols,tMax=tMax,saveAt=saveAt)[1]
            println("Finished Chamber ",i)
        else
            print("Could not complete chamber because previous chamber errored out.")
            sols[i] = ["tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached"]
        end
    end

    return sols

end

"""
    function MultiChamberSolution(tspan::Tuple,u0::Tuple,IV::Vector{Vector{Float64}},IVtimes::Vector{Float64},nChambers::Int64;<keyword arguments>)

Compute the solution for a system of `nChambers` connected CSTRs with `u0` initial conditions. The 
inflow of the first CSTR is given by variable inflow vector `IV` where each entry in `IV` corresponds to the inflow vector at the corresponding entry in `IVtimes`, the outflow of the first CSTR becomes the inflow of second CSTR 
and so on.

This function reads in the parameter values from a .csv file.
The names of the .csv files that contain the parameter values should be "model_parameters.csv", "model_parameters2.csv", 
"model_parameters3.csv", ... and so on.

# Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.
- `saveAt = []`: specific times to save the solution. If given a number `n`, the solver will save the solution every `n` timesteps

# Examples
```julia-repl
julia> u0 = initialConditions();

julia> t = [i for i in 0.0:0.1:50.0];

julia> IV_temp = inflowvector_definition();

julia> IV = [IV_temp*(0.5*rand()+1.0) for i in 1:length(t)];

julia> sols = MultiChamberSolution((0.0,50.0),(u0,u0,u0),IV,t,3);
Finished Chamber 1
Finished Chamber 2
Finished Chamber 3

julia> sols[1]
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 115-element Vector{Float64}:
[...]

u: 115-element Vector{Vector{Float64}}:
[...]

julia> sols[2]
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 115-element Vector{Float64}:
[...]

u: 115-element Vector{Vector{Float64}}:
[...]

julia> sols[3]
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 115-element Vector{Float64}:
[...]

u: 115-element Vector{Vector{Float64}}:
[...]
```
"""
function MultiChamberSolution(tspan::Tuple,u0::Tuple,IV::Vector{Vector{Float64}},IVtimes::Vector{Float64},nChambers::Int64;alg = Rodas4P(), tols=1e-4,tMax = 300.0,saveAt=[])
    
    sols = Vector{Any}(undef,nChambers)

    sols[1] = ADM1sol(tspan,u0[1],IV,IVtimes,alg=alg,tols=tols,tMax=tMax,saveAt=saveAt)[1]

    println("Finished Chamber 1")

    for i in 2:nChambers
        if typeof(sols[i-1]) != Vector{String}
            filename = string("model_parameters",i,".csv")
            sols[i] = ADM1MultiChamberSol(tspan,u0[i],sols[i-1],paramFilename = filename,alg=alg,tols=tols,tMax=tMax,saveAt=saveAt)[1]
            println("Finished Chamber ",i)
        else
            print("Could not complete chamber because previous chamber errored out.")
            sols[i] = ["tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached","tMax Reached"]
        end
    end

    return sols

end

export ExampleMultiChamberSol
"""
    function ExampleMultiChamberSol(tspan::Tuple,u0::Vector,IV::SciMLBase.ODESolution; <keyword arguments>)

Compute the solution for the given timespan, `tspan`; initial condition, `u0`; and inflow vector given by
the solution of the previous CSTR.

Also return the time (in seconds) the solution took to compute.

# Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.

# Examples
```julia-repl
julia> u0 = ADM1code.initialConditions();

julia> IV = ADM1code.inflowvector_definition();

julia> sol, tSol = ADM1code.ExampleMultiChamberSol((0.0,200.0),u0,IV);

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
function ExampleMultiChamberSol(tspan::Tuple,u0::Vector,IV::SciMLBase.ODESolution;alg = Rodas4P(), tols=1e-4,tMax = 300.0)

   # Combine all of the parameter vectors into one vector (parm) for the solver.

   RP = reactorParameterDefinition() # length = 6
   BP = biochemicalparameter_definition() # length = 45
   SP = stoichiometricparameter_definition() # length = 18
   CC = carbonContent_definition() # length = 18
   PhP = physiochemicalParameterDefinition(RP) # length = 20

   parm = zeros(109)
   parm[1:6] = RP
   parm[7:51] = BP
   parm[52:69] = SP
   parm[70:87] = CC
   parm[88:107] = PhP
   parm[108] = time() # start time
   parm[109] = tMax # maximum time

   # the solver algorithm outputs the solution as a vector of solutions at each timestep 
   # convert that to a vector of solutions for each state variable 
   IV_vecs = individualSolutions(IV)

   # convert outflow from previous cstr to gridded interpolation
   global inflowFunctions = [interpolate((IV.t,),IV_vecs[i],Gridded(Constant())) for i in 1:length(IV.u[1])]


   prob = ODEProblem(RHSfunInflowVaried,u0,tspan,parm)

   global sol = "not defined"

   try
      global t = @timed sol = solve(prob,alg, abstol=tols,reltol=tols,isoutofdomain = (u,p,t) -> any(x->x<0,u))
   catch e
      printstyled(stderr,"\nERROR: ", bold=true, color=:red)
      printstyled(stderr,sprint(showerror,e), color=:light_red)
      println(stderr)
   end

   isSolved = (sol != "not defined")

   if isSolved == true
       return sol,t[2]
   else
    #   return nothing # normally the method should return nothing since the problem was not solved
      return ["tMax reached" for i in 1:35],"tMax reached" # only use this for testing
    end

end

export ADM1MultiChamberSol
"""
    function ADM1MultiChamberSol(tspan::Tuple,u0::Vector,IV::SciMLBase.ODESolution; <keyword arguments>)

Compute the solution for the given timespan, `tspan`; initial condition, `u0`; and inflow vector given by
the solution of the previous CSTR.

Also return the time (in seconds) the solution took to compute. The difference between this function and ExampleMultiChamberSol is that this function reads in the parameter values from a .csv file.

# Optional Arguments
- `paramFilename = "model_parameters.csv"`: The name of the .csv that contains the parameter values for this CSTR.
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.
- `saveAt = []`: specific times to save the solution. If given a number `n`, the solver will save the solution every `n` timesteps

# Examples
```julia-repl
julia> u0 = ADM1code.initialConditions();

julia> IV = ADM1code.inflowvector_definition();

julia> sol, tSol = ADM1code.ADM1MultiChamberSol((0.0,200.0),u0,IV);

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
function ADM1MultiChamberSol(tspan::Tuple,u0::Vector,IV::SciMLBase.ODESolution; paramFilename = "model_parameters.csv", alg = Rodas4P(), tols=1e-4,tMax = 300.0,saveAt=[])
   
   # Combine all of the parameter vectors into one vector (parm) for the solver.

   p = CSV.File(paramFilename)
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
         ]


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

   RP = [i[1] for i in RP]
   BP = [i[1] for i in BP]
   SP = [i[1] for i in SP]
   CC = [i[1] for i in CC]
   PhP = [i[1] for i in PhP]

   parm = zeros(109)
   parm[1:6] = RP
   parm[7:51] = BP
   parm[52:69] = SP
   parm[70:87] = CC
   parm[88:107] = PhP
   parm[108] = time() # start time
   parm[109] = tMax # maximum time
   
   # the solver algorithm outputs the solution as a vector of solutions at each timestep 
   # convert that to a vector of solutions for each state variable 
   IV_vecs = individualSolutions(IV)

   # convert outflow from previous cstr to gridded interpolation
   global inflowFunctions = [interpolate((IV.t,),IV_vecs[i],Gridded(Constant())) for i in 1:length(IV.u[1])]

   prob = ODEProblem(RHSfunInflowVaried,u0,tspan,parm)

   global sol = "not defined"

   try
      global t = @timed sol = solve(prob,alg, abstol=tols,reltol=tols,saveat=saveAt,isoutofdomain = (u,p,t) -> any(x->x<0,u))
   catch e
      printstyled(stderr,"\nERROR: ", bold=true, color=:red)
      printstyled(stderr,sprint(showerror,e), color=:light_red)
      println(stderr)
   end

   isSolved = (sol != "not defined")

   if isSolved == true
      return sol,t[2]
   else
      # return nothing # normally the method should return nothing since the problem was not solved
      return ["tMax reached" for i in 1:35],"tMax reached" # only use this for testing
   end

end
