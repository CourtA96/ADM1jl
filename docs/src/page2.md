# Basic Usage

To solve ADM1 for a single tank reactor, just run `ADM1sol`. 

## Getting Started

`ADM1sol` takes the timespan, initial conditions, and inflow vector as inputs. To test this out, run the following code:

```@repl
using ADM1code

u0 = ADM1code.InitialConditions(); # assigns the default initial conditions to u0

IV = ADM1code.inflowvector_definition(); # assigns the default inflow vector to IV

tspan = (0.0,200.0); # the solution will be computed from t=0.0 to t=200.0

sol, tSol = ADM1code.ADM1sol(tspan,u0,IV); # computes the solution  and saves it to sol, the time to solve is saved to tSol

sol # the solution has two fields: t contains the timesteps and u contains the solution at each timestep

tSol # this is the time ExampleSol took to solve the system

```

The initial conditions and inflow vector can be changed to any `Vector{Float64}` of length 35. For a version of the code that is even more flexible, use `ADM1sol`.

## Modifying Parameters

To change the system parameters, such as `T_base` or `P_atm`, just open the file `model_parameters.csv` in the `src` directory. Edit whichever entries are necessary, save, and exit. Running `ADM1sol` again will solve the system with the updated parameters.





