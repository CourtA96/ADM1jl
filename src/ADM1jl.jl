module ADM1jl

#####################################################
# A Julia code to solve the Anaerobic Digestion Model No. 1
# presented by Batstone et al.
#
# Authors:
#  - Alexandra Mazanko (University of Guelph)
#  - Courtney Allen (University of Guelph)
#####################################################

using LinearAlgebra
using DifferentialEquations
using Plots
default(fontfamily="Helvetica")
using SparseArrays
using Memoize
using CSV
using DataFrames
using Interpolations
#####################################################
# Functions to export from other files

export computeDifference
export trueSolutionADM1_200days

export transportmatrix_definition
export petersenmatrixtranspose_definition

export RHSfunInflowVaried
export MultiChamberSolutionExample
export MultiChamberSolution
import ExampleMultiChamberSol
export ADM1MultiChamberSol

export inflowvector_definition
export stoichiometricparameter_definition
export carbonContent_definition
export biochemicalparameter_definition
export reactorParameterDefinition
export physiochemicalParameterDefinition
export computePhysiochemicalParameterDefinition
export initialConditions

export plotSols
export individualSolutions

# dependent files:

include("matrix_definitions.jl")
include("parameter_definitions.jl")
include("plotting_and_splitting.jl")
include("compute_errors.jl")
include("reorder_solutions.jl")
include("compute_steady_state_and_pH.jl")
include("multichamber.jl")

#####################################################
# Functions:

export pressureOfGasses
"""
    pressureOfGasses(sx,php,rp)

Compute the pressures of the gasses.

# Arguments
- `sx::Vector`: the state vector.
- `php::Vector`: the physicochemical parameters.
- `rp::Vector`: the reactor parameters.

# Examples
```jldoctest
julia> u0 = ADM1jl.initialConditions();

julia> rp = ADM1jl.reactorParameterDefinition();

julia> php = ADM1jl.physiochemicalParameterDefinition(rp);

julia> ADM1jl.pressureOfGasses(u0,php,rp)
5-element Vector{Float64}:
    1.6333471490625e-5
    0.6525381992578124
    0.35869584449999997
    1.0669181223042932
    2695.9061152146637
```
"""
function pressureOfGasses(sx,php,rp)

   # PV = nRT or P = SRT where S=n/V := concentration

   p_gas_h2 =  (sx[33] * php[1] * rp[1] / 16)   # S_gas_h2 = sx[33]
   p_gas_ch4 =  (sx[34] * php[1] * rp[1] / 64)  # S_gas_ch4 = sx[34]
   p_gas_co2 =  (sx[35] * php[1] * rp[1])       # S_gas_co2 = sx[35]

   P_gas =  (p_gas_h2 + p_gas_ch4 + p_gas_co2 + php[15])
   q_gas =  (php[16] * (P_gas - rp[3])) # gas flow through pipe eqn (5.10)
   if q_gas < 0
      q_gas = 0
   end

   pressures = [
   p_gas_h2,
   p_gas_ch4,
   p_gas_co2,
   P_gas,
   q_gas
   ]

   return pressures
end

export monod
"""
    monod(u, k)

Compute the monod function `u/(u+k)` where `u` is the state varible and `k` is the half-saturation concentration.

If `u <= 0` return `0`.

# Examples
```jldoctest
julia> ADM1jl.monod(3.0,2.0) # when u is non-zero positive
0.6
```
```jldoctest
julia> ADM1jl.monod(-3.0,2.0) # when u is negative
0
```

"""
function monod(u,k)
   # u is the state variable, k is the half saturation concentration
   if u>0
      return u/(k+u)
   else
      return 0
   end
end

export reactionrates
"""
    reactionrates(bp,rp,php,pressures,sx,NREAC::Int)

Compute and return the vector of reaction rates.

# Arguments
- `bp::Vector`: the biochemical parameters.
- `rp::Vector`: the reactor parameters.
- `php::Vector`: the physiochemical parameters.
- `pressures::Vector`: the gas pressures.
- `sx::Vector`: the state vector.
- `NREAC::Integer`: the number of reaction rates.

# Examples
```jldoctest
julia> u0 = ADM1jl.initialConditions();

julia> bp = ADM1jl.biochemicalparameter_definition();

julia> rp = ADM1jl.reactorParameterDefinition();

julia> php = ADM1jl.physiochemicalParameterDefinition(rp);

julia> pressures = ADM1jl.pressureOfGasses(u0,php,rp);

julia> NREAC = 29;

julia> ADM1jl.reactionrates(bp,rp,php,pressures,u0,NREAC)
29-element Vector{Real}:
 0.155
 0.28
 1.0
 0.29000000000000004
 0.2950855111452082
 ⋮
 0.0
 7.402547081085608e-6
 1.295220255437568
 0.052518812965057435
```
"""
function reactionrates(bp,rp,php,pressures,sx,NREAC::Int)

   "Compute concentration of CO2 and NH4+"

   # S_co2 = S_IC - S_hco3m from eqn (4.3) (ie. S = S_acid+S_base)
   S_co2 = sx[10]-sx[29]

   # S_nh4p = S_IN - S_nh3  since S = S_acid+S_base
   S_nh4p = sx[11]-sx[30]

   "Computing pH"

   K_h20 = php[2]

   # sx[31] = S_cation
   # sx[32] = S_anion

   phi = sx[31]+S_nh4p-sx[29]-(sx[28]/64)-(sx[27]/112)-(sx[26]/160)-(sx[25]/208)-sx[32] # temp variable to solve Eq (B.4) for S_H_pos

   S_H_pos = (-phi + sqrt(phi^2 + 4*K_h20))/2 # solving Eq (B.4) for S_H_pos

   #if typeof(S_H_pos) == Float64
   #   println(S_H_pos,",")
   #end

   """
   For t = (0.0,200.0) the value of S_H_pos starts at 8.321714194473029e-12 and ends at 3.770446836401308e-8. The change is wonky at first, but then begins to increase more slowly.
   """

   #S_H_pos = 0.00000003423

   pH = -log10(S_H_pos)

   "Defining Inhibtion Factors"

   "pH inhibition factors:"

   pHvals = [bp[8] bp[9]; bp[31] bp[32]; bp[37] bp[38]]

   # values of pH_LL and pH_UL
   pH_LLaa = bp[9]
   pH_ULaa = bp[8]
   pH_LLac = bp[32]
   pH_ULac = bp[31]
   pH_LLh2 = bp[38]
   pH_ULh2 = bp[37]

   # parameters
   K_pHaa = 10^(-(pH_LLaa+pH_ULaa)/2.0)
   nn_aa = 3.0/(pH_ULaa-pH_LLaa)
   K_pHac = 10^(-(pH_LLac+pH_ULac)/2.0)
   nn_ac = 3.0/(pH_ULac-pH_LLac)
   K_pHh2 = 10^(-(pH_LLh2+pH_ULh2)/2.0)
   nn_h2 = 3.0/(pH_ULh2-pH_LLh2)

   "pH inhibition factor from BSM2:"
   IpH_aa = (K_pHaa^nn_aa)/(S_H_pos^nn_aa+K_pHaa^nn_aa)
   IpH_ac = (K_pHac^nn_ac)/(S_H_pos^nn_ac+K_pHac^nn_ac)
   IpH_h2 = (K_pHh2^nn_h2)/(S_H_pos^nn_h2+K_pHh2^nn_h2)

   I_pH = [IpH_aa, IpH_ac, IpH_h2]

   "pH inhibition from Table 3.5 c) Empirical upper and lower inhibition"
   #I_pH =(ones(3) + (2 * (10 .^ (0.5*(pHvals[:,2]-pHvals[:,1]))))) ./ (ones(3) + (10 .^ (pH .- pHvals[:,1])) + (10 .^ (pHvals[:,2] .- pH)))

   "pH inhibition from Table 3.5 c) Empirical lower inhibition only"
   #IpH = zeros(3)
   #for i = 1:3
   #  if pH < pHvals[i,1]
   #     IpH[i] = exp(1) .^ (-3*((pH .- pHvals[i,1]) ./ (pHvals[i,1]-pHvals[i,2])) .^ 2)
   #  else
   #     IpH[i] = 1
   #end

   "Non-Competitive Inhibition:"
   # table 3.5 a)
   #KI =  [bp[19], bp[23], bp[27]] # KI = [K_I_h2_fa, K_I_h2_c4, K_I_h2_pro]

   I_h2fa =  (1 / (1 + (sx[8]/ bp[19])))
   I_h2c4 =  (1 / (1 + (sx[8]/ bp[23])))
   I_h2pro =  (1 / (1 + (sx[8]/ bp[27])))

   I_h2 = [I_h2fa,I_h2c4,I_h2pro]

   #I_h2 = 1 ./ (ones(3) + (sx[8] ./ KI))
   if sx[8] < 0
      I_h2 = ones(3)
   end

   I_nh3 = 1 / (1+(sx[30]/bp[33]))

   "Secondary Substrate:"
   # table 3.5 e) where K_S_IN = Ks_NH3_all
   I_IN = (1 / (1+(bp[7]/sx[11])))

   "Final inhibition factors:"

   I = I_pH*I_IN # I_1 in tables 3.1-3.2
   II = I[1]*I_h2 # I_2 in tables 3.1-3.2
   III = I[2]*I_nh3 # I_3 in tables 3.1-3.2

   rrates=Array{Real}(undef,NREAC)
   rrates[1]= bp[1]*sx[13]
   rrates[2]= bp[2]*sx[14]
   rrates[3]= bp[3]*sx[15]
   rrates[4]= bp[4]*sx[16]
   rrates[5]= (bp[10]*(monod(sx[1],bp[11]))*sx[17])*I[1]
   rrates[6]= (bp[13]*(monod(sx[2],bp[14]))*sx[18])*I[1]
   rrates[7]= (bp[16]*(monod(sx[3],bp[17]))*sx[19])*II[1]
   rrates[8]= (bp[20]*(monod(sx[4],bp[21]))*sx[20]*(sx[4]/(sx[4] + sx[5] + 10.0^(-6.0))))*II[2]
   rrates[9]= (bp[20]*(monod(sx[5],bp[21]))*sx[20]*(sx[5]/(sx[5] + sx[4] + 10.0^(-6.0))))*II[2]
   rrates[10]= (bp[24]*(monod(sx[6],bp[25]))*sx[21])*II[3]
   rrates[11]= (bp[28]*(monod(sx[7],bp[29]))*sx[22])*III
   rrates[12]= (bp[34]*(monod(sx[8],bp[35]))*sx[23])*I[3]

   # The decay rates:

   rrates[13] = bp[39]*sx[17]
   rrates[14] = bp[40]*sx[18]
   rrates[15] = bp[41]*sx[19]
   rrates[16] = bp[42]*sx[20]
   rrates[17] = bp[43]*sx[21]
   rrates[18] = bp[44]*sx[22]
   rrates[19] = bp[45]*sx[23]

   "Acid/Base Reactions"

   T_ad = rp[1]
   T_base = rp[2]

   k_AB_va = php[3]
   k_AB_bu = php[4]
   k_AB_pro = php[5]
   k_AB_ac = php[6]
   k_AB_co2 = php[7]
   k_AB_IN = php[8]

   K_ava = php[9]
   K_abu = php[10]
   K_apro = php[11]
   K_aac = php[12]
   K_aco2 = php[13]
   K_aIN = php[14]

   # The rates are given by Equation (5.16)
   # substituting S_base = S-S_acid in rates 20-23
   rrates[20] = k_AB_va*(sx[25]*(S_H_pos+K_ava)-K_ava*sx[4])
   rrates[21] = k_AB_bu*(sx[26]*(S_H_pos+K_abu)-K_abu*sx[5])
   rrates[22] = k_AB_pro*(sx[27]*(S_H_pos+K_apro)-K_apro*sx[6])
   rrates[23] = k_AB_ac*(sx[28]*(S_H_pos+K_aac)-K_aac*sx[7])
   rrates[24] = k_AB_co2*(sx[29]*(S_H_pos+K_aco2)-K_aco2*sx[10])
   rrates[25] = k_AB_IN*(sx[30]*(S_H_pos+K_aIN)-K_aIN*sx[11])

   rrates[26] = 0.0 # cation/Anion

   "Liquid/gas phase transitions"

   p_gas_h2 =  pressures[1]
   p_gas_ch4 =  pressures[2]
   p_gas_co2 =  pressures[3]


   P_gas =  pressures[4]
   q_gas =  pressures[5]

   K_La = php[17]
   K_Hco2 = php[18]
   K_Hch4 = php[19]
   K_Hh2 = php[20]

   rrates[27] = K_La*(sx[8] - 16*K_Hh2*p_gas_h2)
   rrates[28] = K_La*(sx[9] - 64*K_Hch4*p_gas_ch4)
   rrates[29] = K_La*(S_co2 - K_Hco2*p_gas_co2)

   return rrates
end

#function RHSfun(u::Vector{Float64},p::Array{AbstractArray,1},t::Float64)
#function RHSfun(u::Vector{Float64},p::Array{AbstractArray{Float64,N} where N,1},t::Float64)

export RHSfun
"""
    RHSfun(du,u,p,t)

Return the right-hand side of the system of ODEs, this is an in-place function.

# Arguments
- `du::Vector`: the rate change of the state vector (required since the function is defined in-place).
- `u::Vector`: the state vector.
- `p::Vector`: all of the model parameters.
- `t`: the timestep, usually a Float64.
"""
function RHSfun(du,u,p,t)
   # Parameters:

   RP = p[1:6] # reactor parameters
   BP = p[7:51] # biochemical paramters
   SP = p[52:69] # stoichiometric parameters
   CC = p[70:87] # carbon balance parameters
   PhP = p[88:107] # physiochemical parameters
   IV = p[108:142] # inflow vector

   t1 = p[143] # start time
   tMax = p[144] # maximum time

   # compute pressures
   pressures = pressureOfGasses(u,PhP,RP)

   # compute transport matrix
   TM = transportmatrix_definition(RP,pressures)

   # create Petersen Matrix
   PM = petersenmatrixtranspose_definition(RP,BP,SP,CC)

   # compute reaction rates
   rr = reactionrates(BP,RP,PhP,pressures,u,size(PM)[2])

   # compute inflow/outflow rate
   liquidFlow = RP[6]/RP[4]

   t2 = time()

   if t2-t1 > tMax
      errStr = string("Took longer than ",tMax, " seconds.")
      error(errStr)
   end

   du .= TM*u + PM*rr + IV*liquidFlow

end

"""
##########################################################
# MAIN PROGRAM ###########################################
##########################################################
"""

export ExampleSol
"""
    function ExampleSol(tspan::Tuple,u0::Vector,IV::Vector; <keyword arguments>)

Compute the solution for the given timespan, `tspan`; initial condition, `u0`; and inflow vector `IV`.

Also return the time (in seconds) the solution took to compute.

# Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.

# Examples
```julia-repl
julia> u0 = ADM1jl.initialConditions();

julia> IV = ADM1jl.inflowvector_definition();

julia> sol, tSol = ADM1jl.ExampleSol((0.0,200.0),u0,IV);

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
function ExampleSol(tspan::Tuple,u0::Vector,IV::Vector;alg = Rodas4P(), tols=1e-4,tMax = 300.0)

   # Combine all of the parameter vectors into one vector (parm) for the solver.

   RP = reactorParameterDefinition() # length = 6
   BP = biochemicalparameter_definition() # length = 45
   SP = stoichiometricparameter_definition() # length = 18
   CC = carbonContent_definition() # length = 18
   PhP = physiochemicalParameterDefinition(RP) # length = 20
   #IV = inflowvector_definition() # length = 35

   parm = zeros(144)
   parm[1:6] = RP
   parm[7:51] = BP
   parm[52:69] = SP
   parm[70:87] = CC
   parm[88:107] = PhP
   parm[108:142] = IV
   parm[143] = time() # start time
   parm[144] = tMax # maximum time

   prob=ODEProblem(RHSfun,u0,tspan,parm)

   global sol = "not defined"

   try
      global t = @timed sol = solve(prob,alg, abstol=tols,reltol=tols)
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

export ADM1sol
"""
    function ADM1sol(tspan::Tuple,u0::Vector,IV::Vector; <keyword arguments>)

Compute the solution for the given timespan, `tspan`; initial condition, `u0`; and inflow vector `IV`.

Also return the time (in seconds) the solution took to compute. The difference between this function and ExampleSol is that this function reads in the parameter values from a .csv file.

# Optional Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.

# Examples
```julia-repl
julia> u0 = initialConditions();

julia> IV = inflowvector_definition();

julia> sol, tSol = ADM1sol((0.0,200.0),u0,IV);

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
function ADM1sol(tspan::Tuple,u0::Vector,IV::Vector{Float64};alg = Rodas4P(), tols=1e-4,tMax = 300.0)
   
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

   RP = [i[1] for i in RP]
   BP = [i[1] for i in BP]
   SP = [i[1] for i in SP]
   CC = [i[1] for i in CC]
   PhP = [i[1] for i in PhP]

   parm = zeros(144)
   parm[1:6] = RP
   parm[7:51] = BP
   parm[52:69] = SP
   parm[70:87] = CC
   parm[88:107] = PhP
   parm[108:142] = IV
   parm[143] = time() # start time
   parm[144] = tMax # maximum time

   prob=ODEProblem(RHSfun,u0,tspan,parm)

   global sol = "not defined"

   try
      global t = @timed sol = solve(prob,alg, abstol=tols,reltol=tols)
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

"""
    function ADM1sol(tspan::Tuple,u0::Vector,IV::Vector{Vector{Float64}},IVtimes::Vector{Float64}; <keyword arguments>)

Compute the solution with variable inflow for the given timespan, `tspan`; initial condition, `u0`; and variable inflow vector 
`IV` where each entry in IV corresponds to the inflow vector at the corresponding entry in `IVtimes`

Also return the time (in seconds) the solution took to compute.

# Optional Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.

# Examples
```julia-repl
julia> u0 = ADM1jl.initialConditions();

julia> t = [i for i in 0.0:0.1:50.0]

julia> IV_temp = ADM1jl.inflowvector_definition();

julia> IV = [IV_temp*(0.5*rand()+1.0) for i in 1:length(t)]

julia> sol,tSol = ADM1jl.ADM1sol((0.0,50.0),u0,IV,t);

julia> sol
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 2011-element Vector{Float64}:
[...]

u: 2011-element Vector{Vector{Float64}}:
[...]

julia> tSol
14.4643811

```
"""
function ADM1sol(tspan::Tuple,u0::Vector,IV::Vector{Vector{Float64}},IVtimes::Vector{Float64};alg = Rodas4P(), tols=1e-4,tMax = 300.0)
   
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

   parm = zeros(144)
   parm[1:6] = RP
   parm[7:51] = BP
   parm[52:69] = SP
   parm[70:87] = CC
   parm[88:107] = PhP
   parm[108] = time() # start time
   parm[109] = tMax # maximum time

   IV_vecs = [[IV[i][j] for i in 1:length(IV)] for j in 1:35]

   global inflowFunctions = [interpolate((IVtimes,),IV_vecs[i],Gridded(Constant())) for i in 1:length(IV[1])]

   prob=ODEProblem(RHSfunInflowVaried,u0,tspan,parm)

   global sol = "not defined"

   try
      global t = @timed sol = solve(prob,alg, abstol=tols,reltol=tols)
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

end # module
