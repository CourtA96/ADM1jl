# matrix_definitions.jl

"""
The functions that build the:
   - transport matrix
   - petersen matrix
   - inflow vector
"""

export initialConditions_biofilm
function initialConditions_biofilm()
   """
   Function that returns the same initial conditions as used in PyADM1.py
   """

   # Order of variables:
   # [
   # S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I,
   # X_xc, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I,
   # S_va_ion, S_bu_ion, S_pro_ion, S_ac_ion, S_hco3_ion, S_nh3, S_cat, S_an
   # S_gas_h2, S_gas_ch4, S_gas_co2
   # ]

   u0 = [
   0.012, 0.0053, 0.099, 0.012, 0.013, 0.016, 0.2, 2.3e-7, 0.055, 0.15, 0.13, 0.33, # substrates
   0.31, 0.028, 0.1, 0.029, 0.42, 1.18, 0.24, 0.43, 0.14, 0.76, 0.32, 25.6, # suspended biomass
   0.011, 0.013, 0.016, 0.2, 0.14, 0.0041, 0.040, 0.020, # ions
   1.02e-5, 1.63, 0.014, # gas
   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 # biofilm biomass
   ]

   return u0

end

export transportmatrix_definition_biofilm
@memoize function transportmatrix_definition_biofilm(rp,pressures)
   """
   #####################################################
   # this is where the input of the transport matrix takes place
   # ---------------------------------------------------
   # returns:
   #   TM: transport matrix sparse array with dilution rates
   #####################################################
   """

   liquidFlow = rp[6]/rp[4]
   gasFlow = pressures[5]/rp[5]

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

   TM_temp = [TMtemp;zeros(12)]

   TM = sparse(1:47,1:47,-TM_temp)

   return TM
end

export petersenmatrixtranspose_definition_biofilm
@memoize function petersenmatrixtranspose_definition_biofilm(rp,bp,sp,cc)
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
   # P2 contains the column of the Petersen Matrix (Rows of Transpose)
   P2 = [
   1, 1, 1, 1, 1, 1,
   2, 2, 2, 2,
   3, 3, 3, 3,
   4, 4, 4, 4,
   5, 5, 5, 5, 5, 5,
   6, 6, 6, 6, 6, 6, 6, 6,
   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
   8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
   9, 9, 9, 9, 9,
   10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
   11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
   12, 12,
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
   35,
   36, 36, 36, 36, 36, 36, 36, 36, # Biofilm
   37, 37,
   38, 38, 
   39, 39, 
   40, 40, 
   41, 41, 
   42, 42,
   43, 43, 43,
   44, 44,
   45, 45,
   46, 46,
   47
   ]

   # P1 contains the row of the Petersen Matrix (Columns of Transpose)
   P1 = [
   2, 4, 5, 31, 33, 34,  #1
   3, 6, 32, 35, #2
   4, 7, 33, 36, #3
   6, 8, 35, 37, #4
   5, 6, 9, 34, 35, 38, #5
   5, 6, 8, 10, 34, 35, 37, 39, #6
   5, 6, 7, 8, 9, 10, 11, 34, 35, 36, 37, 38, 39, 40, #7
   5, 6, 7, 8, 9, 10, 12, 27, 34, 35, 36, 37, 38, 39, 41, #8
   11, 12, 28, 40, 41, #9
   5, 6, 10, 11, 12, 29, 34, 35, 39, 40, 41,  #10
   1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, #11 (Rates 1, 13-19 are from BSM2)
   30, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, #11 (continued)
   1, 30, #12
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
   29, #35
   30, 42, 43, 44, 45, 46, 47, 48, #36 Biofilm
   30, 31, #37
   30, 32, #38
   30, 33, #39
   34, 42, #40
   35, 43, #41
   36, 44, #42
   37, 38, 45, #43
   39, 46, #44
   40, 47, #45
   41, 48, #46
   30, #47
   ]

   # Q contains the entries in the Petersen Matrix
   Q = [
   1.0, 1.0-f_fali, -1.0, 1.0, 1.0-f_fali, -1.0, #1
   1.0, -1.0, 1.0, -1.0, #2
   f_fali, -1.0, f_fali, -1.0, #3
   (1.0-Y_aa)*f_vaaa, -1.0, (1.0-Y_aa)*f_vaaa, -1.0, #4
   (1.0-Y_su)*f_busu, (1.0-Y_aa)*f_buaa, -1.0, (1.0-Y_su)*f_busu, (1.0-Y_aa)*f_buaa, -1.0, #5
   (1.0-Y_su)*f_prosu, (1.0-Y_aa)*f_proaa, (1.0-Y_c4)*0.54, -1.0, (1.0-Y_su)*f_prosu, (1.0-Y_aa)*f_proaa, (1.0-Y_c4)*0.54, -1.0, #6
   (1.0-Y_su)*f_acsu, (1.0-Y_aa)*f_acaa, (1.0-Y_fa)*0.7, (1.0-Y_c4)*0.31, (1.0-Y_c4)*0.8, (1.0-Y_pro)*0.57, -1.0, (1.0-Y_su)*f_acsu, (1.0-Y_aa)*f_acaa, (1.0-Y_fa)*0.7, (1.0-Y_c4)*0.31, (1.0-Y_c4)*0.8, (1.0-Y_pro)*0.57, -1.0, #7
   (1.0-Y_su)*f_h2su, (1.0-Y_aa)*f_h2aa, (1.0-Y_fa)*0.3, (1.0-Y_c4)*0.15, (1.0-Y_c4)*0.2, (1.0-Y_pro)*0.43, -1.0, -1.0, (1.0-Y_su)*f_h2su, (1.0-Y_aa)*f_h2aa, (1.0-Y_fa)*0.3, (1.0-Y_c4)*0.15, (1.0-Y_c4)*0.2, (1.0-Y_pro)*0.43, -1.0, #8
   (1.0-Y_ac), (1.0-Y_h2), -1.0, (1.0-Y_ac), (1.0-Y_h2), #9
   -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, #10
   (N_xc- f_xIxc*N_I - f_sIxc*N_I - f_prxc*N_aa), -Y_su*Nbac, N_aa-(Y_aa*Nbac), -Y_fa*Nbac, -Y_c4*Nbac, -Y_c4*Nbac, -Y_pro*Nbac, -Y_ac*Nbac, -Y_h2*Nbac, #11 (part 1)
   (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), #11 (part 2)
   (N_xc- f_xIxc*N_I - f_sIxc*N_I - f_prxc*N_aa), -Y_su*Nbac, N_aa-(Y_aa*Nbac), -Y_fa*Nbac, -Y_c4*Nbac, -Y_c4*Nbac, -Y_pro*Nbac, -Y_ac*Nbac, -Y_h2*Nbac, #11 (part 3)
   (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), (Nbac-N_xc), #11 (part 4)
   f_sIxc, f_sIxc, #12 from ADM1
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
   V_liq/V_gas, #35
   -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, #36 Biofilm
   f_chxc, -1.0, #37
   f_prxc, -1.0, #38
   f_lixc, -1.0, #39
   Y_su, -1.0, #40
   Y_aa, -1.0, #41
   Y_fa, -1.0, #42
   Y_c4, Y_c4, -1.0, #43
   Y_pro, -1.0, #44
   Y_ac, -1.0, #45
   Y_h2, -1.0, #46
   f_xIxc, #47 from ADM1
   ]

   # Q = convert(Array{Real},Q)

   PM = sparse(P2,P1,Q) # matrix transpose

   """Carbon balances"""
   C = zeros(P2[end])
   # Ctemp = zeros(P2[end])
   # C = convert(Array{Real},Ctemp)
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
   C[36:47] = C[13:24]

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
   PM[10,30] = -PM[:,30]'*C # Biofilm
   PM[10,31] = -PM[:,31]'*C
   PM[10,32] = -PM[:,32]'*C
   PM[10,33] = -PM[:,33]'*C
   PM[10,34] = -PM[:,34]'*C
   PM[10,35] = -PM[:,35]'*C
   PM[10,36] = -PM[:,36]'*C
   PM[10,37] = -PM[:,37]'*C
   PM[10,38] = -PM[:,38]'*C
   PM[10,39] = -PM[:,39]'*C
   PM[10,40] = -PM[:,40]'*C
   PM[10,41] = -PM[:,41]'*C
   PM[10,42] = -PM[:,42]'*C
   PM[10,43] = -PM[:,43]'*C
   PM[10,44] = -PM[:,44]'*C
   PM[10,45] = -PM[:,45]'*C
   PM[10,46] = -PM[:,46]'*C
   PM[10,47] = -PM[:,47]'*C
   PM[10,48] = -PM[:,48]'*C

   return Matrix(PM)
end

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

export reactionrates_biofilm
"""
    reactionrates_biofilm(bp,rp,php,pressures,sx,NREAC::Int)

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
julia> u0 = ADM1code.InitialConditions();

julia> bp = ADM1code.biochemicalparameter_definition();

julia> rp = ADM1code.reactorParameterDefinition();

julia> php = ADM1code.physiochemicalParameterDefinition(rp);

julia> pressures = ADM1code.pressureOfGasses(u0,php,rp);

julia> NREAC = 48;

julia> ADM1code.reactionrates_biofilm(bp,rp,php,pressures,u0,NREAC)
48-element Vector{Real}:
[...]
```
"""
function reactionrates_biofilm(bp,rp,php,pressures,sx,NREAC::Int)

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


   "Biofilm Rates"
   
   rrates[30]= bp[1]*sx[36]
   rrates[31]= bp[2]*sx[37]
   rrates[32]= bp[3]*sx[38]
   rrates[33]= bp[4]*sx[39]
   rrates[34]= (bp[10]*(monod(sx[1],bp[11]))*sx[40])*I[1]
   rrates[35]= (bp[13]*(monod(sx[2],bp[14]))*sx[41])*I[1]
   rrates[36]= (bp[16]*(monod(sx[3],bp[17]))*sx[42])*II[1]
   rrates[37]= (bp[20]*(monod(sx[4],bp[21]))*sx[43]*(sx[4]/(sx[4] + sx[5] + 10.0^(-6.0))))*II[2]
   rrates[38]= (bp[20]*(monod(sx[5],bp[21]))*sx[43]*(sx[5]/(sx[5] + sx[4] + 10.0^(-6.0))))*II[2]
   rrates[39]= (bp[24]*(monod(sx[6],bp[25]))*sx[44])*II[3]
   rrates[40]= (bp[28]*(monod(sx[7],bp[29]))*sx[45])*III
   rrates[41]= (bp[34]*(monod(sx[8],bp[35]))*sx[46])*I[3]

   # The decay rates:

   rrates[42] = bp[39]*sx[40]
   rrates[43] = bp[40]*sx[41]
   rrates[44] = bp[41]*sx[42]
   rrates[45] = bp[42]*sx[43]
   rrates[46] = bp[43]*sx[44]
   rrates[47] = bp[44]*sx[45]
   rrates[48] = bp[45]*sx[46]

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
   TM = transportmatrix_definition_biofilm(RP,pressures)

   # create Petersen Matrix
   PM_temp = petersenmatrixtranspose_definition_biofilm(RP,BP,SP,CC)
   PM_bio = petersenmatrixtranspose_biofilmAdditions()

   PM = [PM_temp PM_bio]

   # compute reaction rates
   rr_temp = reactionrates_biofilm(BP,RP,PhP,pressures,u,48)
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

Compute the biofilm solution for the given timespan, `tspan`; initial condition, `u0`; and inflow vector `IV`.
The last 12 elements of the statevector are 

Also return the time (in seconds) the solution took to compute.

# Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.

# Examples
```julia-repl
julia> u0 = initialConditions_biofilm(); # 47 elements

julia> IV = inflowvector_definition();

julia> sol, tSol = Biofilm((0.0,200.0),u0,IV);

julia> sol
retcode: Success
Interpolation: specialized 2nd order "free" stiffness-aware interpolation
t: 778-element Vector{Float64}:
[...]

u: 778-element Vector{Vector{Float64}}:
[...]

julia> tSol
11.4968794
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

   try
      global t = @timed sol = solve(prob,alg=Rosenbrock23())
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
      return ["tMax reached" for i in 1:47],"tMax reached" # only use this for testing
   end

end

export Biofilm
"""
    function Biofilm(tspan::Tuple,u0::Vector,IV::Vector{Vector{Float64}},IVtimes::Vector{Float64}; <keyword arguments>)

Compute the biofilm solution for the given timespan, `tspan`; initial condition, `u0`; and variable inflow vector `IV`.
The last 12 elements of the statevector are 

Also return the time (in seconds) the solution took to compute.

# Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.

# Examples
```julia-repl
julia> u0 = initialConditions_biofilm(); # 47 elements

julia> IV = inflowvector_definition();

julia> sol, tSol = Biofilm((0.0,200.0),u0,IV);

julia> sol
retcode: Success
Interpolation: specialized 2nd order "free" stiffness-aware interpolation
t: 778-element Vector{Float64}:
[...]

u: 778-element Vector{Vector{Float64}}:
[...]

julia> tSol
11.4968794
```
"""
function Biofilm(tspan::Tuple,u0::Vector,IV::Vector{Vector{Float64}},IVtimes::Vector{Float64}; tols=1e-4, tMax = 300.0, saveAt=[])

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

   IV_vecs = [[IV[i][j] for i in 1:length(IV)] for j in 1:35]

   global inflowFunctions = [interpolate((IVtimes,),IV_vecs[i],Gridded(Linear())) for i in 1:length(IV[1])]

   prob=ODEProblem(biofilmRHSfunInflowVaried,u0,tspan,parm)

   global sol = "not defined"

   try
      global t = @timed sol = solve(prob,alg=Rosenbrock23())
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
      return ["tMax reached" for i in 1:47],"tMax reached" # only use this for testing
   end

end

export biofilmRHSfunInflowVaried
"""
    biofilmRHSfunInflowVaried(du,u,p,t)

Return the right-hand side of the system of 0D biofilm ODEs when inflow is varied, this is an in-place function.

# Arguments
- `du::Vector`: the rate change of the state vector (required since the function is defined in-place).
- `u::Vector`: the state vector.
- `p::Vector`: all of the model parameters.
- `t`: the timestep, usually a Float64.
"""
function biofilmRHSfunInflowVaried(du,u,p,t)
   # Parameters:

   RP = p[1:6] # reactor parameters
   BP = p[7:51] # biochemical paramters
   SP = p[52:69] # stoichiometric parameters
   CC = p[70:87] # carbon balance parameters
   PhP = p[88:107] # physiochemical parameters
   IV_temp = [inflowFunctions[i][t] for i in 1:35] # inflow is the value of the interpolation of the previous cstr at time t

   IV = [IV_temp; zeros(12)]

   t1 = p[143] # start time
   tMax = p[144] # maximum time

   AD = p[145:146] # attach and detach params

   # compute pressures
   pressures = pressureOfGasses(u,PhP,RP)

   # compute transport matrix
   TM = transportmatrix_definition_biofilm(RP,pressures)

   # create Petersen Matrix
   PM_temp = petersenmatrixtranspose_definition_biofilm(RP,BP,SP,CC)
   PM_bio = petersenmatrixtranspose_biofilmAdditions()

   PM = [PM_temp PM_bio]

   # compute reaction rates
   rr_temp = reactionrates_biofilm(BP,RP,PhP,pressures,u,48)
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

   du .= TM*u + PM*rr - TM*IV

end

export ADM1toBiofilmSolution
"""
    function ADM1toBiofilmSolution(tspan::Tuple,u0::Tuple,IV::Vector,nChambers::Int64; <keyword arguments>)

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
function ADM1toBiofilmSolution(tspan::Tuple,u0::Tuple,IV::Vector;alg = Rodas4P(), tols=1e-4,tMax = 300.0,saveAt=[])
    
    sols = Vector{Any}(undef,2)

    sols[1] = ADM1sol(tspan,u0[1],IV,alg=alg,tols=tols,tMax=tMax,saveAt=saveAt)[1]

    println("Finished Chamber 1")

    filename = string("model_parameters2.csv")
    sols[2] = BiofilmMultiChamberSol(tspan,u0[2],sols[1],paramFilename = filename,tols=tols,tMax=tMax,saveAt=saveAt)[1]
    println("Finished Chamber 2")

    return sols

end

export MultiBiofilmSolution
"""
    function MultiBiofilmSolution(tspan::Tuple,u0::Tuple,IV::Vector,nChambers::Int64; <keyword arguments>)

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
function MultiBiofilmSolution(tspan::Tuple,u0::Tuple,IV::Vector,nChambers::Int64;alg = Rodas4P(), tols=1e-4,tMax = 300.0,saveAt=[])
    
    sols = Vector{Any}(undef,nChambers)

    sols[1] = BiofilmMultiChamberSol(tspan,u0[1],IV,tols=tols,tMax=tMax,saveAt=saveAt)[1]

    println("Finished Chamber 1")

    for i in 2:nChambers
        filename = string("model_parameters",i,".csv")
        sols[i] = BiofilmMultiChamberSol(tspan,u0[i],sols[i-1],paramFilename = filename,tols=tols,tMax=tMax,saveAt=saveAt)[1]
        println("Finished Chamber ",i)
    end

    return sols

end

export BiofilmMultiChamberSol
"""
    function BiofilmMultiChamberSol(tspan::Tuple,u0::Vector,IV::Vector{Vector{Float64}},IVtimes::Vector{Float64}; <keyword arguments>)

Compute the biofilm solution for the given timespan, `tspan`; initial condition, `u0`; and inflow vector `IV`.
The last 12 elements of the statevector are 

Also return the time (in seconds) the solution took to compute.

# Arguments
- `alg = Rodas4P()`: the ODE solver algorithm.
- `tols = 1e-4`: the absolute and relative tolerance of the solver method.
- `tMax = 300.0`: the maximum time (in seconds), that the function will run before timing out.

# Examples
```julia-repl
julia> u0 = initialConditions_biofilm(); # 47 elements

julia> IV = inflowvector_definition();

julia> sol, tSol = Biofilm((0.0,200.0),u0,IV);

julia> sol
retcode: Success
Interpolation: specialized 2nd order "free" stiffness-aware interpolation
t: 778-element Vector{Float64}:
[...]

u: 778-element Vector{Vector{Float64}}:
[...]

julia> tSol
11.4968794
```
"""
function BiofilmMultiChamberSol(tspan::Tuple,u0::Vector,IV::SciMLBase.ODESolution; paramFilename = "model_parameters.csv", tols=1e-4, tMax = 300.0, saveAt=[])

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
   # parm[108:142] = IV
   parm[143] = time() # start time
   parm[144] = tMax # maximum time
   parm[145:146] = AD

   IV_vecs = individualSolutions(IV)

   global inflowFunctions = [interpolate((IV.t,),IV_vecs[i],Gridded(Constant())) for i in 1:length(IV.u[1])]

   prob=ODEProblem(biofilmRHSfunInflowVaried,u0,tspan,parm)

   global sol = "not defined"

   try
      global t = @timed sol = solve(prob,alg=Rosenbrock23())
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
      return ["tMax reached" for i in 1:47],"tMax reached" # only use this for testing
   end

end

export plotSolsBiofilm
function plotSolsBiofilm(sol;titleText::String="Plots of Solutions")
   # I got the code to create the title from https://stackoverflow.com/questions/43066957/adding-global-title-to-plots-jl-subplots
   y = (ones(3))
   title1 = scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], text(string(titleText, " 1"))),axis=false, grid=false, leg=false,size=(200,100),reuse=false)
   title2 = scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], text(string(titleText, " 2"))),axis=false, grid=false, leg=false,size=(200,100),reuse=false)
   title3 = scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], text(string(titleText, " 3"))),axis=false, grid=false, leg=false,size=(200,100),reuse=false)


   sols = individualSolutions(sol)
   pltNums = [i for i in 1:20]
   figure1temp = plot([(sol.t,sols[i]) for i in 1:20],layout=20,title = ["($i)" for j in 1:1, i in pltNums], titleloc = :right,reuse=false, titlefont = font(8),labels=:none)
   figure1 = plot(title1,figure1temp,layout=grid(2,1,heights=[0.05,0.9]),size=(1200,800))
   display(figure1)

   pltNums = [i for i in 21:40]
   figure2temp = plot([(sol.t,sols[i]) for i in 21:40],layout=20,title = ["($i)" for j in 1:1, i in pltNums], titleloc = :right,reuse=false, titlefont = font(8),labels=:none)
   figure2 = plot(title2,figure2temp,layout=grid(2,1,heights=[0.05,0.9]),size=(1200,800))
   display(figure2)

   pltNums = [i for i in 41:size(sols)[1]]
   figure3temp = plot([(sol.t,sols[i]) for i in 41:47],layout=20,title = ["($i)" for j in 1:1, i in pltNums], titleloc = :right,reuse=false, titlefont = font(8),labels=:none)
   figure3 = plot(title3,figure3temp,layout=grid(2,1,heights=[0.05,0.9],size=(1200,800)))
   display(figure3)
end