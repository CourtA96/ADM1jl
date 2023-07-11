#commpute_steady_state_and_pH.jl

function computepH(soln)
   pHvec = []

   R =  0.083145
   T_base =  298.15
   T_ad =  308.15

   K_h2o = (1e-14) * exp((55900 / (100 * R)) * (1 / T_base - 1 / T_ad))

   for i in 1:length(soln.t)
      sol = soln[i]
      # S_co2 = S_IC - S_hco3m from eqn (4.3) (ie. S = S_acid+S_base)
      S_co2 = sol[10]-sol[29]
      # S_nh4p = S_IN - S_nh3  since S = S_acid+S_base
      S_nh4p = sol[11]-sol[30]

      phi = sol[31]+S_nh4p-sol[29]-(sol[28]/64)-(sol[27]/112)-(sol[26]/160)-(sol[25]/208)-sol[32] # temp variable to solve Eq (B.4) for S_H_pos
      S_H_pos = (-phi + sqrt(phi^2 + 4*K_h2o))/2

      pH = -log10(S_H_pos)
      pHvec = vcat(pHvec,pH)
   end
   return pHvec
end

function compute_du(u,IV)
   ######################################################
   # defines right hand side for chemostat models:
   # ----------------------------------------------------
   # note:  normally the frst argument should be
   #        du, but this seems to not work here, but if
   #        it is omitted, it runs
   # ----------------------------------------------------
   # on input:
   #   u:  state variable
   #   p:  array of arrays containg all model parameters
   #       p[1]: transpose of Petersen matrix
   #       p[2]: inflow vector
   #       p[3]: vector of biochemical parameters
   #       p[4]: vector of reactor parameters
   #       p[5]: vector of physiochemical parameters
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

   RP = reactorParameterDefinition() # length = 6
   BP = biochemicalparameter_definition() # length = 45
   SP = stoichiometricparameter_definition() # length = 18
   CC = carbonContent_definition() # length = 18
   PhP = physiochemicalParameterDefinition(RP) # length = 15

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

   du = TM*u + PM*rr + IV*liquidFlow

end

function timeToSteadyState(sol;IV=ADM1code.inflowvector_definition(),tol=1e-4)
   N = length(sol)
   v = zeros(N-1)

   for i in 1:N-1
      v[i] = abs(norm(compute_du(sol[i],IV))/(sol.t[i+1]-sol.t[i]))
   end

   ind = findall(x->x<=tol,v)

   if length(ind)>0
      indexSteady = minimum(ind)
      return sol.t[indexSteady]
   else
      return "Does not converge"
   end

end

function computeChange(sol;IV=ADM1code.inflowvector_definition())

   N = length(sol)
   v = zeros(N-1)

   for i in 1:N-1
      v[i] = abs(norm(compute_du(sol[i],IV))/(sol.t[i+1]-sol.t[i]))
   end

   return (v[1],v[end])

end
