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

Transport matrix

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

   TM = sparse(1:35,1:35,-TMtemp)

   return TM
end

export petersenmatrixtranspose_definition
"""
    petersenmatrixtranspose_definition(rp,bp,sp,cc)

Transpose of the Petersen matrix

```
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
   # P2 contains the column of the Petersen Matrix
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

   # P1 contains the row of the Petersen Matrix
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

   # Q contains the entries in the Petersen Matrix
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

   # Q = convert(Array{Real},Q)

   PM = sparse(P2,P1,Q)

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
