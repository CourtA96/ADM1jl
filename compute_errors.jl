#compute_errors.jl

function PyADM1results()
   """
   Function that returns a vector containing the solutions given by PyADM1.py at time t=200.0
   """
   results = [
   0.01199361878146199, 0.005200381337872014, 0.09898093692189579, 0.011868933687433882, 0.013074922196440932, 0.01597024556830165, 0.199862345903217, 2.3399607909219572e-07, 0.05530940017982466, 0.15001903585035922, 0.12999908025333287, 0.32999999810480457, #substrates
   0.3099927282048765, 0.027996040159098214, 0.10025714796109751, 0.02904994126445399, 0.4200010468463959, 1.179987011031539, 0.24000062802816663, 0.43000064239933733, 0.14000101865957973, 0.7600029896614757, 0.3199983928491184, 25.600010410706354, # biomass
   0.011839659664365938, 0.013045392220216352, 0.015928914690047917, 0.1994695337234892, 0.14029205509292164, 0.004084592387200918, 0.04, 0.02, # ions
   1.022845431128427e-05, 1.6352974486422807, 0.013906264769753896 # gas
   ]

   return results

end

function computeError(sol)
   """
   Function that returns a vector containing the relative difference between the input sol and the result of PyADM1.py.

   Note:
   The PyADM1.py solution is at the time t=200.0
   """
   trueSol = PyADM1results()
   error = abs.((sol[end]-trueSol)./trueSol)*100
   return error
end

function computeError2(sol,trueSol)
   """
   Function that returns a vector containing the relative difference between the inputs sol and trueSol
   """
   error = abs.((sol[end]-trueSol)./trueSol)*100
   return error
end
