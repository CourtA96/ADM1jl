# Getting started

To solve ADM1 for a single tank reactor, just run `ADM1sol`. For multiple tanks running in parallel, use `multichamberSolution`

## Installation

To install `ADM1jl` open the Julia REPL and run the following:

```
julia> using Pkg
julia> Pkg.add("https://github.com/CourtA96/ADM1jl")
```

OR enter the Pkg REPL by pressing ] and running:

```
add https://github.com/CourtA96/ADM1jl
```

Then, to use the package, run

```
julia> using ADM1jl
```


## Getting Started

Before beginning, make sure that the file `model_parameters.csv` is saved in your working directory. `model_parameters.csv` can be found on the github [here](https://github.com/CourtA96/ADM1jl/blob/main/model_parameters.csv). 

`ADM1sol` takes the timespan, initial conditions, and inflow vector as inputs. The timespan is length 2 and type `Tuple{Float64}`. It specifies how the initial and final times of the simulation. The initial conditions and inflow vector both have type `Vector{Float64}` and length 35.  To test this out, run the following code:

```@repl
using ADM1jl

u0 = ADM1jl.InitialConditions(); # assigns the default initial conditions to u0

IV = ADM1jl.inflowvector_definition(); # assigns the default inflow vector to IV

tspan = (0.0,200.0); # the solution will be computed from t=0.0 to t=200.0

sol, tSol = ADM1jl.ADM1sol(tspan,u0,IV); # computes the solution  and saves it to sol, the time to solve is saved to tSol

sol # the solution has two fields: t contains the timesteps and u contains the solution at each timestep

tSol # this is the time ExampleSol took to solve the system

```

The initial conditions and inflow vector can be changed to any `Vector{Float64}` of length 35. For a version of the code that is even more flexible, use `ADM1sol`.

## Modifying Parameters

To change the system parameters, such as `T_base` or `P_atm`, just open the file `model_parameters.csv` in the `src` directory. Edit whichever entries are necessary, save, and exit. Running `ADM1sol` again will solve the system with the updated parameters.

## State Variables and their Indices

The follow table shows which state variable corresponds to each index.

| Index | State Variable |
| :---: | :------------: |
| 1     | `S_su`         |
| 2     | `S_aa`         |
| 3     | `S_fa`         |
| 4     | `S_va`         |
| 5     | `S_bu`         |
| 6     | `S_pr`         |
| 7     | `S_ac`         |
| 8     | `S_h2`         |
| 9     | `S_ch4`        |
| 10    | `S_IC`         |
| 11    | `S_IN`         |
| 12    | `S_I`          |
| 13    | `X_xc`         |
| 14    | `X_ch`         |
| 15    | `X_pr`         |
| 16    | `X_li`         |
| 17    | `X_su`         |
| 18    | `X_aa`         |
| 19    | `X_fa`         |
| 20    | `X_c4`         |
| 21    | `X_pro`        |
| 22    | `X_ac`         |
| 23    | `X_h2`         |
| 24    | `X_I`          |
| 25    | `S_va_ion`     |
| 26    | `S_bu_ion`     |
| 27    | `S_pro_ion`    |
| 28    | `S_ac_ion`     |
| 29    | `S_hco3_ion`   |
| 30    | `S_nh3`        |
| 31    | `S_cat`        |
| 32    | `S_an`         |
| 33    | `S_gas_h2`     |
| 34    | `S_gas_ch4`    |
| 35    | `S_gas_co2`    |