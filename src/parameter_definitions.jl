#parameter_definitions.jl

"""
The functions that define the model parameters:
   - inflow vector
   - stoichiometic parameters
   - carbon content
   - biochemical parameters
   - reactor parameters
   - physiochemical parameters
   - initial conditions
"""

"""
      inflowvector_definition()

Returns the inflow vector.

# Examples
```jldoctest
julia> ADM1code.inflowvector_definition()
35-element Vector{Float64}:
 0.01
 0.001
 0.001
 0.001
 ⋮
 0.02
 0.0
 0.0
 0.0
```
"""
@memoize function inflowvector_definition()
   """
   #####################################################
   # this is where the input of the inflow vector takes place
   # this follows a simple array type
   # on input:
   #     IV: inflow vector information
   # ---------------------------------------------------
   # returns:
   #   IV: single column array containing the inflow vector
   #####################################################
   """

   # Order of variables:
   S_su = 0.01; S_aa = 0.001; S_fa = 0.001; S_va = 0.001; S_bu = 0.001;
   S_pro = 0.001; S_ac = 0.001; S_h2 = 1e-8; S_ch4 = 1e-5; S_IC = 0.04;
   S_IN = 0.01; S_I = 0.02;

   X_xc = 2.0; X_ch = 5.0; X_pr = 20.0; X_li = 5.0; X_su = 0.0; X_aa = 0.01;
   X_fa = 0.01; X_c4 = 0.01; X_pro = 0.01; X_ac = 0.01; X_h2 = 0.01; X_I = 25.0;

   S_va_ion = 0.0; S_bu_ion = 0.0; S_pro_ion = 0.0; S_ac_ion = 0.0; S_hco3_ion = 0.0; S_nh3 = 0.0; S_cat = 0.04; S_an = 0.02;

   S_gas_h2 = 0.0; S_gas_ch4 = 0.0; S_gas_co2 = 0.0;

   IV = [
      S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I,
      X_xc, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I,
      S_va_ion, S_bu_ion, S_pro_ion, S_ac_ion, S_hco3_ion, S_nh3, S_cat, S_an,
      S_gas_h2, S_gas_ch4, S_gas_co2
      ]

   return IV
end

@memoize function stoichiometricparameter_definition()
   """
      In the order that they are presented in Table 6.1 of ADM1
   """
   f_slxc = 0.1; f_xlxc = 0.20; f_chxc = 0.20; f_prxc = 0.20; f_lixc = 0.3; #f_lixc = 0.25;
   N_xc = 0.002; N_I = 0.002; f_fali = 0.95;
   f_h2su = 0.19; f_busu = 0.13; f_prosu = 0.27; f_acsu = 0.41; f_h2aa = 0.06;
   N_aa = 0.007; f_vaaa = 0.23; f_buaa = 0.26; f_proaa = 0.05; f_acaa = 0.40

   SP =  [
   f_slxc, f_xlxc, f_chxc, f_prxc, f_lixc,
   N_xc, N_I, f_fali,
   f_h2su, f_busu, f_prosu, f_acsu, f_h2aa,
   N_aa, f_vaaa, f_buaa, f_proaa, f_acaa
   ]

   return SP
end

@memoize function carbonContent_definition()
   """
      In order of i. Values are given in Table 2.6 of ADM1
   """
   C_su =  0.0313; C_aa =  0.0300; C_fa =  0.0217; C_va =  0.0240;
   C_bu =  0.025; C_prop =  0.0268; C_ac =  0.0313; C_h2 = 0.0; C_ch4 = 0.0156;
   C_IC = 1.0; C_IN = 0.0; C_xc = 0.02786; C_ch = 0.0313; C_protein = 0.03; C_li = 0.022;
   C_sI = 0.03; C_xI = 0.03; C_bac = 0.0313

   CC = [
   C_su, C_aa, C_fa, C_va, C_bu, C_prop, C_ac, C_h2, C_ch4, C_IC, C_IN, C_sI,
   C_xc, C_ch, C_protein, C_li, C_bac, C_xI
   ]
   return CC
end

@memoize function biochemicalparameter_definition()
   """
    In the order that they are presented in Table 6.2 of ADM1
   """
    kdis = 0.5; khyd_CH = 10.0; khyd_PR = 10.0; khyd_LI = 10.0; tresX = 40.0;
    kdec_all = 0.02; Ks_NH3_all = 0.0001;
    pHUL_acet_acid = 5.5; pHLL_acet_acid = 4.0;
    km_su = 30.0; KS_SU = 0.5; Ysu = 0.10;
    km_aa = 50.0; KS_aa = 0.3; Yaa =0.08;
    km_fa = 6.0; KS_fa = 0.4; Yfa = 0.06; KlH2_fa = 5.0e-6;
    km_c4 = 20.0; KS_c4 = 0.2; Yc4 = 0.06; KlH2_c4 = 1e-5;
    km_pro = 13.0; KS_pro = 0.1; Ypro = 0.04; KlH2_pro = 3.5e-6;
    km_ac = 8.0; KS_ac = 0.15; Yac = 0.05; pHUL_ac = 7.0; pHLL_ac = 6.0; KINH3 = 0.0018;
    km_h2 = 35.0; KS_h2 = 7.0e-6;
    Yh2 = 0.06; pHUL_h2 = 6.0; pHLL_h2 = 5.0;
    k_decSU = 0.02; k_decAA = 0.02; k_decFA = 0.02; k_decC4 = 0.02; k_decPRO = 0.02; k_decAC = 0.02; k_decH2 = 0.02;

   BP = [kdis, khyd_CH, khyd_PR, khyd_LI, tresX,
         kdec_all, Ks_NH3_all,
         pHUL_acet_acid, pHLL_acet_acid,
         km_su, KS_SU, Ysu,
         km_aa, KS_aa, Yaa,
         km_fa, KS_fa, Yfa, KlH2_fa,
         km_c4, KS_c4, Yc4, KlH2_c4,
         km_pro, KS_pro, Ypro, KlH2_pro,
         km_ac, KS_ac, Yac, pHUL_ac, pHLL_ac, KINH3,
         km_h2, KS_h2,
         Yh2, pHUL_h2, pHLL_h2,
         k_decSU, k_decAA, k_decFA, k_decC4, k_decPRO, k_decAC, k_decH2
         ]

   return BP
end

@memoize function reactorParameterDefinition()
   """
   List of parameters that relate to the reactor itself.
   From
   """
   # Temperature
   T_ad = 308.15 # Temperature of reactor Potential issue here for BSM2Defaults
   T_base = 298.15 # Base temperature for parameters. Should not be changed.

   # Atmospheric pressure
   P_atm = 1.013

   # Volumes of liquid and gas
   V_liq =  3400
   V_gas =  300

   # Flow Rate
   Q_ad = 170.0

   RP = [
   T_ad, T_base,
   P_atm,
   V_liq, V_gas,
   Q_ad
   ]

   return RP

end

@memoize function physiochemicalParameterDefinition(rp::Vector{Float64})
   """
   Parameters from tables 4.1 and 4.2. Some parameters are not specified by ADM1 and are taken from the BSM2, which can be found at: http://iwa-mia.org/benchmarking/#BSM2

   Equation 4.10 is used to calculate temperature variations.
   """
   R = 0.083145 # L*bar*K^(-1)*mol^(-1)

   # In the following equations R is multiplied by 100 to convert to SI units (J*K^(-1)*mol^(-1))

   k_AB_va = 1.0e10
   k_AB_bu = 1.0e10
   k_AB_pro = 1.0e10
   k_AB_ac = 1.0e10
   k_AB_co2 = 1.0e10
   k_AB_IN = 1.0e10

   # From Table 4.1
   K_h2o = (1e-14)*exp(55900.0*((1/rp[2]) - (1/rp[1]))/(100.0*R)) # Eq (4.10)
   K_ava = 10^(-4.86)
   K_abu = 10^(-4.82)
   K_apro = 10^(-4.88)
   K_aac = 10^(-4.76)
   K_aco2 = 10^(-6.35)*exp(7646.0*((1/rp[2]) - (1/rp[1]))/(100.0*R)) # Eq (4.10)
   K_aIN = 10^(-9.25)*exp(51965.0*((1/rp[2]) - (1/rp[1]))/(100.0*R)) # Eq (4.10)

   p_gas_h2o = 0.0313 * exp(5290.0 * ((1.0/rp[2]) - (1.0/rp[1]))) # Eq (5.7)
   kp = 5.0e4

   K_La = 200.0 # from BSM2 (pg 38 of ADM1)

   # Table 4.2
   K_Hco2 = 0.035*exp(-19410.0*((1/rp[2]) - (1.0/rp[1]))/(100.0*R)) # Eq (4.10)
   K_Hch4 = 0.0014*exp(-14240.0*((1/rp[2]) - (1.0/rp[1]))/(100.0*R)) # Eq (4.10)
   K_Hh2 = 7.8e-4*exp(-4180.0*((1/rp[2]) - (1.0/rp[1]))/(100.0*R)) # Eq (4.10)

   PhP = [
   R,
   K_h2o, k_AB_va, k_AB_bu, k_AB_pro, k_AB_ac, k_AB_co2, k_AB_IN,
   K_ava, K_abu, K_apro, K_aac, K_aco2, K_aIN,
   p_gas_h2o, kp,
   K_La, K_Hco2, K_Hch4, K_Hh2,
   ]

   return PhP
end

@memoize function computePhysiochemicalParameterDefinition(rp::Vector{Float64},php::Vector{Float64})
   """
   Overloaded method.

   Calculates PhP using inputs:
      - rp:    reactor parameters
      - php:   independent physiochemical parameters
   Parameters from tables 4.1 and 4.2. Some parameters are not specified by ADM1 and are taken from the BSM2, which can be found at: http://iwa-mia.org/benchmarking/#BSM2

   Equation 4.10 is used to calculate temperature variations.
   """
   R = php[1] # L*bar*K^(-1)*mol^(-1)

   # In the following equations R is multiplied by 100 to convert to SI units (J*K^(-1)*mol^(-1))

   k_AB_va = 10^php[4]
   k_AB_bu = 10^php[5]
   k_AB_pro = 10^php[6]
   k_AB_ac = 10^php[7]
   k_AB_co2 = 10^php[8]
   k_AB_IN = 10^php[9]

   # From Table 4.1
   K_h2o = 10^php[2]*exp(php[3]*((1/rp[2]) - (1/rp[1]))/(100.0*R)) # Eq (4.10)
   K_ava = 10^php[10]
   K_abu = 10^php[11]
   K_apro = 10^php[12]
   K_aac = 10^php[13]
   K_aco2 = 10^php[14]*exp(php[15]*((1/rp[2]) - (1/rp[1]))/(100.0*R)) # Eq (4.10)
   K_aIN = 10^php[16]*exp(php[17]*((1/rp[2]) - (1/rp[1]))/(100.0*R)) # Eq (4.10)

   p_gas_h2o = php[18] * exp(php[19] * ((1.0/rp[2]) - (1.0/rp[1]))) # Eq (5.7)
   kp = php[20]

   K_La = php[21] # from BSM2 (pg 38 of ADM1)

   # Table 4.2
   K_Hco2 = php[22]*exp(php[23]*((1/rp[2]) - (1.0/rp[1]))/(100.0*R)) # Eq (4.10)
   K_Hch4 = php[24]*exp(php[25]*((1/rp[2]) - (1.0/rp[1]))/(100.0*R)) # Eq (4.10)
   K_Hh2 = php[26]*exp(php[27]*((1/rp[2]) - (1.0/rp[1]))/(100.0*R)) # Eq (4.10)

   PhP = [
   R,
   K_h2o, k_AB_va, k_AB_bu, k_AB_pro, k_AB_ac, k_AB_co2, k_AB_IN,
   K_ava, K_abu, K_apro, K_aac, K_aco2, K_aIN,
   p_gas_h2o, kp,
   K_La, K_Hco2, K_Hch4, K_Hh2,
   ]

   return PhP
end

function InitialConditions()
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
   0.31, 0.028, 0.1, 0.029, 0.42, 1.18, 0.24, 0.43, 0.14, 0.76, 0.32, 25.6, # biomass
   0.011, 0.013, 0.016, 0.2, 0.14, 0.0041, 0.040, 0.020, # ions
   1.02e-5, 1.63, 0.014 # gas
   ]

   return u0

end
