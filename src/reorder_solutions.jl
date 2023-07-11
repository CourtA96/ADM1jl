#reorder_solutions.jl

function printVec(vec::Vector)
    print("[",vec[1],)
    for v in vec[2:end]
       print(",\n",v)
    end
    print("]")
 end

function reorderFromPython(sol)
   return vcat(sol[1:24],sol[28:32],sol[34],sol[25:26],sol[36:38])
end

function reorderFromJava(sol)
   return vcat(sol[1:24],sol[27:32],sol[25:26],sol[33:35])
end

function reorderToPython(sol)
   # S_co2 = S_IC - S_hco3m from eqn (4.3) (ie. S = S_acid+S_base)
   S_co2 = sol[10]-sol[29]
   # S_nh4p = S_IN - S_nh3  since S = S_acid+S_base
   S_nh4p = sol[11]-sol[30]

   R =  0.083145
   T_base =  298.15
   T_ad =  308.15

   K_h2o = (1e-14) * exp((55900 / (100 * R)) * (1 / T_base - 1 / T_ad))

   phi = sol[31]+S_nh4p-sol[29]-(sol[28]/64)-(sol[27]/112)-(sol[26]/160)-(sol[25]/208)-sol[32] # temp variable to solve Eq (B.4) for S_H_pos
   S_H_pos = (-phi + sqrt(phi^2 + 4*K_h2o))/2

   return vcat(sol[1:24],sol[31:32],S_H_pos,sol[25:29],S_co2,sol[30],S_nh4p,sol[33:35])
end

function reorderToJava(sol)
   # S_co2 = S_IC - S_hco3m from eqn (4.3) (ie. S = S_acid+S_base)
   S_co2 = sol[10]-sol[29]
   # S_nh4p = S_IN - S_nh3  since S = S_acid+S_base
   S_nh4p = sol[11]-sol[30]

   R =  0.083145
   T_base =  298.15
   T_ad =  308.15

   K_h2o = (1e-14) * exp((55900 / (100 * R)) * (1 / T_base - 1 / T_ad))

   phi = sol[31]+S_nh4p-sol[29]-(sol[28]/64)-(sol[27]/112)-(sol[26]/160)-(sol[25]/208)-sol[32] # temp variable to solve Eq (B.4) for S_H_pos
   S_H_pos = (-phi + sqrt(phi^2 + 4*K_h2o))/2

   pH = -log10(S_H_pos)

   return vcat(sol[1:24],sol[31,32],sol[25:30],sol[33:35],[170.0,35.0,0.0,0.0],pH,S_co2,S_nh4p)
end

function reorderInflowToPython(IV)
   return vcat(IV[1:24],IV[31:32])
end
