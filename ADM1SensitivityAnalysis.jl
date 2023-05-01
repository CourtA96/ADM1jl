module ADM1SensitivityAnalysis

function RHSfun2(du,u,p,t)
   ######################################################
   # defines right hand side for chemostat models:
   # ----------------------------------------------------
   # note:  normally the frst argument should be
   #        du, but this seems to not work here, but if
   #        it is omitted, it runs
   # ----------------------------------------------------
   # on input:
   #   u:  state variable
   #   p:  array of arrays containing all model parameters to be changed
   #        - A 1D array of Int64 that contains the indices of the paramters
   #          must be globally defined outside of the function.
   #   t:  time (not used when system is autonomous)
   # returns:
   #   du:  rate of change for state variable
   ######################################################
   # N.B.:  i determine the number of process here from
   #    the size of the petersen matrix as size(p[2])[2].
   #    explicitly handing the number of processes over
   #    to the function reactionrates made programming it
   #    much easier
   ######################################################

   # Parameters:

   parms = zeros(126)
   parms[1:6] = rpd
   parms[7:51] = bpd
   parms[52:68] = spd
   parms[69:76] = ccd
   parms[77:91] = ppd
   parms[92:126] = ivd

   parms[parmsToVary] = p

   RP = parms[1:6] # reactor parameters
   BP = parms[7:51] # biochemical paramters
   SP = parms[52:68] # stoichiometric parameters
   CC = parms[69:76] # carbon balance parameters
   PhP = parms[77:91] # physiochemical parameters
   IV = parms[92:126] # inflow vector

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

   du .= convert(Array{Real},TM*u + PM*rr + IV*liquidFlow)

end

function RHSfun3(du,u,p,t)
   ######################################################
   # defines right hand side for chemostat models:
   # ----------------------------------------------------
   # note:  normally the frst argument should be
   #        du, but this seems to not work here, but if
   #        it is omitted, it runs
   # ----------------------------------------------------
   # on input:
   #   u:  state variable
   #   p:  A 2D array of Float64 containing the parameters to be varied and
   #       their indices
   #        - p[1,:] contains the paramters
   #        - p[2,:] contains their indices (which are converted to Int64
   #          inside the function)
   #   t:  time (not used when system is autonomous)
   # returns:
   #   du:  rate of change for state variable
   ######################################################
   # N.B.:  i determine the number of process here from
   #    the size of the petersen matrix as size(p[2])[2].
   #    explicitly handing the number of processes over
   #    to the function reactionrates made programming it
   #    much easier
   ######################################################

   # Parameters:

   parms = zeros(126)
   parms[1:6] = reactorParameterDefinition()
   parms[7:51] = biochemicalparameter_definition()
   parms[52:68] = stoichiometricparameter_definition()
   parms[69:76] = carbonContent_defintion()
   parms[77:91] = physiochemicalParameterDefinition(parms[1:6])
   parms[92:126] = inflowvector_definition()

   # convert indices from Float64 to Int64 and set specified entries of parm to be equal to p[1,:]
   parms[Int.(p[2,:])] = p[1,:]

   RP = parms[1:6] # reactor parameters
   BP = parms[7:51] # biochemical paramters
   SP = parms[52:68] # stoichiometric parameters
   CC = parms[69:76] # carbon balance parameters
   PhP = parms[77:91] # physiochemical parameters
   IV = parms[92:126] # inflow vector

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

   du .= convert(Array{Real},TM*u + PM*rr + IV*liquidFlow)

end

function LocalSensitivitySol(tspan::Tuple,u0::Vector)
    """
    # definition of model parameters
    #-------------------------------
    # for now parameters are defined ad hoc, eventually
    # they should be provided externally
    #---------------------------------------------------
    # TM: transport matrix
    # PM: transpose of Petersen matrix
    # iv: inflow vector
    # rp: vector of reaction parameters
    #
    # parm:  combines the above into one for ODEsolver
    """

    RP = reactorParameterDefinition() # length = 6
    BP = biochemicalparameter_definition() # length = 45
    SP = stoichiometricparameter_definition() # length = 17
    CC = carbonContent_defintion() # length = 8
    PhP = physiochemicalParameterDefinition(RP) # length = 15
    IV = inflowvector_definition() # length = 35

    parm = zeros(126)
    parm[1:6] = RP
    parm[7:51] = BP
    parm[52:68] = SP
    parm[69:76] = CC
    parm[77:91] = PhP
    parm[92:126] = IV

    prob = ODEForwardSensitivityProblem(RHSfun, u0, tspan, parm)

    alg = Rodas5()
    sol = solve(prob,alg)

    return sol
 end

 function RHSfun2Test(tspan::Tuple,u0::Vector)
    """
    # Takes a Tuple specifing the desired timespan and a Vector containing the
    # initial values.
    #
    # Computes the solution using RHSfun2
    # parm:        contains the parameters that are varied
    # parmsToVary: contains the indices of the parameters to be varied
    """

    parm = zeros(3)
    parm[1] = 1.02
    parm[2] = 3400
    parm[3] = 400

    global parmsToVary = [3,4,5]

    global rpd = reactorParameterDefinition()
    global bpd = biochemicalparameter_definition()
    global spd = stoichiometricparameter_definition()
    global ccd = carbonContent_defintion()
    global ppd = physiochemicalParameterDefinition(rpd)
    global ivd = inflowvector_definition()

    prob = ODEProblem(RHSfun2,u0,tspan,parm)

    alg = Rosenbrock23()
    sol = solve(prob,alg)

    return sol
 end

 function RHSfun3Test(tspan::Tuple,u0::Vector)
    """
    # Takes a Tuple specifing the desired timespan and a Vector containing the
    # initial values.
    # parm:  contains the parameters that are varied and their indices
    #  - parm[1,:] contains the parameters
    #  - parm[2,:] contains the indices (these are input as Float64, and are
    #    converted to Int64 inside RHSfun3)
    """

    parm = zeros(2,3)
    parm[:,1] = [1.02, 3.0]
    parm[:,2] = [3400.0, 4.0]
    parm[:,3] = [400.0, 5.0]

    prob = ODEProblem(RHSfun3,u0,tspan,parm)

    alg = Rosenbrock23()
    sol = solve(prob,alg)

    return sol
 end

end #module
