# matrix_definitions.jl

"""
The functions that build the:
   - transport matrix
   - petersen matrix
   - inflow vector
"""

export transportmatrix_definition
"""
    transportmatrix_definition(rp,pressures)

The transport matrix describes the flow of material out of the system.

# Arguments
- `rp::Vector`: reactor parameters
- `pressures::Vector`: pressures computed by pressureOfGasses()

```
"""
@memoize function transportmatrix_definition(rp,pressures)
   """
   #####################################################
   # this is where the input of the transport matrix takes place
   # ---------------------------------------------------
   # returns:
   #   TM: transport matrix sparse array with dilution rates
   #####################################################
   """

   liquidFlow = rp[6]/rp[4] # liquid flow rate (Q_ad/V_liq)
   gasFlow = pressures[5]/rp[5] # gas flow rate (q_gas/V_gas) (q_gas depends on k_p the pipe resistence coefficient)

   TMtemp = [
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   liquidFlow,
   0.0,
   0.0,
   0.0,
   0.0,
   0.0,
   0.0,
   liquidFlow,
   liquidFlow,
   gasFlow,
   gasFlow,
   gasFlow
   ]

   TM = sparse(1:35,1:35,-TMtemp)

   return TM
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
```juia-repl
julia> u0 = initialConditions();

julia> bp = biochemicalparameter_definition();

julia> rp = reactorParameterDefinition();

julia> php = physiochemicalParameterDefinition(rp);

julia> pressures = pressureOfGasses(u0,php,rp);

julia> NREAC = 29;

julia> reactionrates(bp,rp,php,pressures,u0,NREAC)
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

   # If the number of reaction rates is changed, the number of 
   # rows in the Petersen matrix transpose must be changed accordingly

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

   "pH inhibition factors from ADM1:" # replaced with BSM2 inhibition factors
   #"pH inhibition from Table 3.5 c) Empirical upper and lower inhibition"
   #I_pH =(ones(3) + (2 * (10 .^ (0.5*(pHvals[:,2]-pHvals[:,1]))))) ./ (ones(3) + (10 .^ (pH .- pHvals[:,1])) + (10 .^ (pHvals[:,2] .- pH)))

   #"pH inhibition from Table 3.5 c) Empirical lower inhibition only"
   #IpH = zeros(3)
   #for i = 1:3
   #  if pH < pHvals[i,1]
   #     IpH[i] = exp(1) .^ (-3*((pH .- pHvals[i,1]) ./ (pHvals[i,1]-pHvals[i,2])) .^ 2)
   #  else
   #     IpH[i] = 1
   #end

   "Non-Competitive Inhibition:"
   # table 3.5 a) of ADM1
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

   I = I_pH*I_IN # I_1 in tables 3.1-3.2 of ADM1
   II = I[1]*I_h2 # I_2 in tables 3.1-3.2 of ADM1
   III = I[2]*I_nh3 # I_3 in tables 3.1-3.2 of ADM1

   "Reaction Rates"

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

   # The rates are given by Equation (5.16) in ADM1
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

export petersenmatrixtranspose_definition
"""
    petersenmatrixtranspose_definition(rp,bp,sp,cc)

Transpose of the Petersen matrix

# Arguments
- `rp::Vector`: reactor parameters
- `bp::Vector`: biochemical parameters
- `sp::Vector`: stoichiometric parameters
- `cc::Vector`: carbon content

"""
@memoize function petersenmatrixtranspose_definition(rp,bp,sp,cc)
   """
   #####################################################
   # this is where the input of the peteresen matrix takes place
   # this follows a coordinate system where
   # on input:
   #     RP : reactor parameters
   #     P1: row coordinates for paramters
   #     P2: column coordintes for parameters
   #     Q: parameter values for corresponding coordinates
   # ---------------------------------------------------
   # returns:
   #   PM: transposed petersen matrix of a sparse array form
   #####################################################
   """

   """Liquid/gas phase transitions"""
   # Parameters from http://iwa-mia.org/benchmarking/#BSM2
   V_liq =  rp[4]
   V_gas =  rp[5]

   """ Stoichiometric Parameters """

   f_sIxc = sp[1]; f_xIxc = sp[2]; f_chxc = sp[3]; f_prxc = sp[4]; f_lixc = sp[5];
   N_xc = sp[6]; N_I = sp[7]; f_fali = sp[8];
   f_h2su = sp[9]; f_busu = sp[10]; f_prosu = sp[11]; f_acsu = sp[12]; f_h2aa = sp[13];
   N_aa = sp[14]; f_vaaa = sp[15]; f_buaa = sp[16]; f_proaa = sp[17]; f_acaa = sp[18]

   """ Other Parameters """
   Y_su = bp[12]
   Y_aa = bp[15]
   Y_fa = bp[18]
   Y_c4 = bp[22]
   Y_pro = bp[26]
   Y_ac = bp[30]
   Y_h2 = bp[36]
   Nbac = 0.08/14

   """Constructing Petersen matrix"""
   # P2 contains the columns of the Petersen Matrix
   P2 = [
   1, 1, 1,
   2, 2,
   3, 3,
   4, 4,
   5, 5, 5,
   6, 6, 6, 6,
   7, 7, 7, 7, 7, 7, 7,
   8, 8, 8, 8, 8, 8, 8, 8,
   9, 9, 9,
   10, 10, 10, 10, 10, 10,
   11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
   12,
   13, 13, 13, 13, 13, 13, 13, 13,
   14, 14,
   15, 15,
   16, 16,
   17, 17,
   18, 18,
   19, 19,
   20, 20, 20,
   21, 21,
   22, 22,
   23, 23,
   24,
   25, # acids
   26,
   27,
   28,
   29,
   30,
   31, # cation/anion
   32,
   33, # gas phase
   34,
   35
   ]

   # P1 contains the rows of the Petersen Matrix (corresponding column in commented on right)
   P1 = [
   2, 4, 5, #1
   3, 6, #2
   4, 7, #3
   6, 8, #4
   5, 6, 9, #5
   5, 6, 8, 10, #6
   5, 6, 7, 8, 9, 10, 11, #7
   5, 6, 7, 8, 9, 10, 12, 27, #8
   11, 12, 28, #9
   5, 6, 10, 11, 12, 29, #10
   1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, #11 (Rates 1, 13-19 are from BSM2)
   1, #12
   1, 13, 14, 15, 16, 17, 18, 19, #13
   1, 2, #14
   1, 3, #15
   1, 4, #16
   5, 13, #17
   6, 14, #18
   7, 15, #19
   8, 9, 16, #20
   10, 17, #21
   11, 18, #22
   12, 19, #23
   1, #24
   20, #25 (acids)
   21, #26
   22, #27
   23, #28
   24, #29
   25, #30
   26, #31 (cation/anion)
   26, #32
   27, #33 gas phase equations
   28, #34
   29 #35
   ]

   # Q contains the entries of the Petersen Matrix
   Q = [
   1.0, 1.0-f_fali, -1.0, #1
   1.0, -1.0, #2
   f_fali, -1.0, #3
   (1.0-Y_aa)*f_vaaa, -1.0, #4
   (1.0-Y_su)*f_busu, (1.0-Y_aa)*f_buaa, -1.0, #5
   (1.0-Y_su)*f_prosu, (1.0-Y_aa)*f_proaa, (1.0-Y_c4)*0.54, -1.0, #6
   (1.0-Y_su)*f_acsu, (1.0-Y_aa)*f_acaa, (1.0-Y_fa)*0.7, (1.0-Y_c4)*0.31, (1.0-Y_c4)*0.8, (1.0-Y_pro)*0.57, -1.0, #7
   (1.0-Y_su)*f_h2su, (1.0-Y_aa)*f_h2aa, (1.0-Y_fa)*0.3, (1.0-Y_c4)*0.15, (1.0-Y_c4)*0.2, (1.0-Y_pro)*0.43, -1.0, -1.0, #8
   (1.0-Y_ac), (1.0-Y_h2), -1.0, #9
   -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, #10
   (N_xc- f_xIxc*N_I - f_sIxc*N_I - f_prxc*N_aa), -Y_su*Nbac, N_aa-(Y_aa*Nbac), -Y_fa*Nbac, -Y_c4*Nbac, -Y_c4*Nbac, -Y_pro*Nbac, -Y_ac*Nbac, -Y_h2*Nbac, #11 (part 1)
   (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), #11 (part 2)
   f_sIxc, #12 from ADM1
   -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, #13
   f_chxc, -1.0, #14
   f_prxc, -1.0, #15
   f_lixc, -1.0, #16
   Y_su, -1.0, #17
   Y_aa, -1.0, #18
   Y_fa, -1.0, #19
   Y_c4, Y_c4, -1.0, #20
   Y_pro, -1.0, #21
   Y_ac, -1.0, #22
   Y_h2, -1.0, #23
   f_xIxc, #24 from ADM1
   -1.0, #25  (acids)
   -1.0, #26
   -1.0, #27
   -1.0, #28
   -1.0, #29
   -1.0, #30
   1.0, #31 (cation/anion)
   -1.0, #32
   V_liq/V_gas, #33 (gas phase equations)
   V_liq/V_gas, #34
   V_liq/V_gas #35
   ]

   PM = sparse(P2,P1,Q)

   """Carbon balances"""
   C = zeros(P2[end])

   C[1] =  cc[1] #SU kmole C.kg^-1COD
   C[2] =  cc[2] #aa kmole C.kg^-1COD
   C[3] =  cc[3] #fa kmole C.kg^-1COD
   C[4] =  cc[4] #va kmole C.kg^-1COD
   C[5] =  cc[5] #bu kmole C.kg^-1COD
   C[6] =  cc[6] #pro kmole C.kg^-1COD
   C[7] =  cc[7] #ac kmole C.kg^-1COD
                 #h2
   C[9] =  cc[9] #ch4 kmole C.kg^-1COD
                 #IC
                 #IN
   C[12] = cc[12] #sI
   C[13] = cc[13] #xc
   C[14] = cc[14] #ch
   C[15] = cc[15] #protein
   C[16] = cc[16] #li
   C[17:23] = cc[17]*ones(7) #bac
   C[24] = cc[18] # aa

   PM[10,1] = -PM[:,1]'*C
   PM[10,2] = -PM[:,2]'*C
   PM[10,3] = -PM[:,3]'*C
   PM[10,4] = -PM[:,4]'*C
   PM[10,5] = -PM[:,5]'*C
   PM[10,6] = -PM[:,6]'*C
   PM[10,7] = -PM[:,7]'*C
   PM[10,8] = -PM[:,8]'*C
   PM[10,9] = -PM[:,9]'*C
   PM[10,10] = -PM[:,10]'*C
   PM[10,11] = -PM[:,11]'*C
   PM[10,12] = -PM[:,12]'*C
   PM[10,13] = -PM[:,13]'*C
   PM[10,14] = -PM[:,14]'*C
   PM[10,15] = -PM[:,15]'*C
   PM[10,16] = -PM[:,16]'*C
   PM[10,17] = -PM[:,17]'*C
   PM[10,18] = -PM[:,18]'*C
   PM[10,19] = -PM[:,19]'*C

   return Matrix(PM)
end
