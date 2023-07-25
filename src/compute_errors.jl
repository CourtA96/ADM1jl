#compute_errors.jl
"""
The Functions that 
   - define the true solution at 200 days 
   - compute the difference between two input solutions
"""

export computeDifference
"""
    computeDifference(sol, trueSol)

Compute the percent difference between each element in two vector of type `Vector{Float64}`.

# Arguments
- `sol::Vector`: the solution given at a time.
- `trueSol::Vector`: the true solution at that time.

# Examples
```jldoctest
julia> u0 = initialConditions();

julia> IV = inflowvector_definition();

julia> sol, tSol = ADM1sol((0.0,200.0),u0,IV); # compute the solution

julia> trueSol = trueSolutionADM1_200days; # the true solution at 200 days

julia> d = computeDifference(sol[end],trueSol); # compute the difference
```
"""
function computeDifference(sol,trueSol)
   """
   Function that returns a vector containing the relative difference between the inputs sol and trueSol
   """
   error = abs.((sol-trueSol)./trueSol)*100
   return error
end

export trueSolutionADM1_200days
"""
    trueSolutionADM1_200days()

Return the solution (`Vector{Float64}`) for the default parameters at 200 days.
```
"""
function trueSolutionADM1_200days()
   return [0.011954829810563401, 0.0053147401832121036, 0.09862141249882697, 0.011625006887037806, 0.013250730223207484, 0.015783663823179507, 0.19738633671962316, 2.3594504387169125e-7, 0.05508920869990237, 0.15268433594303998, 0.1302294760799933, 0.3286981032420448, 0.30869798826346184, 0.027947243664334303, 0.10257410933597609, 0.029483054551227825, 0.420166031332133, 1.1791718405506475, 0.24303536337252007, 0.431921107655618, 0.13730593572985716, 0.7605713976442028, 0.3170229998944603, 25.61739467100561, 0.011596224203413465, 0.013220740490142689, 0.015742812872396397, 0.1969985369743189, 0.1427840760901829, 0.004087783314163287, 0.04, 0.02, 1.0241037663368444e-5, 1.6256194037384155, 0.014150346745679485]
end
