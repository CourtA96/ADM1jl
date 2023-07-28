# Getting started

To solve ADM1 for a single tank reactor, just run `ADM1sol`. For multiple tanks running in parallel, use `multichamberSolution`


## Installation

To install `ADM1jl` open the Julia REPL and run the following:

```
julia> using Pkg
julia> Pkg.add(url="https://github.com/CourtA96/ADM1jl")
```

OR enter the Pkg REPL by pressing ] and running:


```
add https://github.com/CourtA96/ADM1jl
```

Then, to use the package, run

```
julia> using ADM1jl
```


## Basic Usage

Before beginning, make sure that the file `model_parameters.csv` is saved in your working directory. `model_parameters.csv` can be found on the github [here](https://github.com/CourtA96/ADM1jl/blob/main/model_parameters.csv). 

`ADM1sol` takes the timespan, initial conditions, and inflow vector as inputs. The timespan is length 2 and type `Tuple{Float64}`. It specifies how the initial and final times of the simulation. The initial conditions and inflow vector both have type `Vector{Float64}` and length 35. To test this out, run the following code:

```@repl
using ADM1jl

u0 = initialConditions(); # assigns the default initial conditions to u0

IV = inflowvector_definition(); # assigns the default inflow vector to IV

tspan = (0.0,200.0); # the solution will be computed from t=0.0 to t=200.0

sol, tSol = ADM1sol(tspan,u0,IV); # computes the solution  and saves it to sol, the time to solve is saved to tSol

sol # the solution has two fields: t contains the timesteps and u contains the solution at each timestep

tSol # this is the time ExampleSol took to solve the system

```

### Modifying Parameters

To change the system parameters, such as `T_base` or `P_atm`, just open the file `model_parameters.csv` in your working directory (`model_parameters.csv` can be found [here](https://github.com/CourtA96/ADM1jl/blob/main/model_parameters.csv)). Edit whichever entries are necessary, save, and exit. Running `ADM1sol` again will solve the system with the updated parameters.

### Specifying Alorithms

By default, `ADM1sol` solves the system using the `Rodas4P` algorithm given in the `DifferentialEquations` package (documentation for `DifferentialEquations` available [here](https://docs.sciml.ai/DiffEqDocs/stable/)). To use a different algorithm, install the `DifferentialEquations` package and follow the example below:

```@repl
using DifferentialEquations

using ADM1jl

u0 = initialConditions();

IV = inflowvector_definition();

sol,tSol = ADM1sol((0.0,50.0),u0,IV, alg = Rosenbrock23()); # solve the system using the Rosenbrock23 algorithm
```
In principle, any ODE solver listed in the `DifferentialEquations` [documentation](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) can be used to solve the system. However, ADM1 is a stiff system of equations, so many solvers may not be stable.

### Variable Inflow

To solve the system for variable inflow, specify the inflow as a `Vector{Vector{Float}}`, ie. as a vector that contains the vector of inflow conditions at a different times, these times are indexed by an additional input vector `t`.

The `ADM1sol` function interpolates the function using the `interpolate` function from the `Interpolations` package (documentation for `Interpolations` available [here](https://juliamath.github.io/Interpolations.jl/stable/control/)). The `Gridded(Linear())` interpolation algorithm is specified, meaning that the interpolation between timesteps is linear. Interpolation allows for the use of adaptive step methods, which are more stable than fixed stepsize methods. As when the inflow is fixed, the default solver algorithm is `Rodas4P`.

The following code breaks the timespan into 0.1 day increments and then randomizes the inflow concentrations within 50 percent of the default values every 0.1 days. The system is then solved for these randomly varied inflow vectors.

```@repl
using ADM1jl

u0 = initialConditions(); # the default initial conditions

tspan = (0.0,50.0); # The timespan of the solution is 50 days

t = [i for i in tspan[1]:0.1:tspan[2]]; # Break the timespan into 0.1 day increments

IV_temp = inflowvector_definition(); # temporary inflow vector

IV = [IV_temp*(0.5*rand()+1.0) for i in 1:length(t)] # for each time in t, vary the inflow concentrations within 50 percent of their default values

sol,tSol = ADM1sol((0.0,50.0),u0,IV,t); # solve the system

sol # The solution

tSol # the time it took to solve

```

## Multiple Reactors in Parallel

To model multiple reactors in series, (ie. where the outflow from the first reactor becomes to inflow to the second, and so on) make sure that there are `model_parameters.csv` files corresponding to each reactor in your working directory. These files should be called `model_parameters.csv`, `model_parameters2.csv`, `model_parameters3.csv`, and so on.

To solve multiple reactors, use the `MultiChamberSolution` function, which takes the timespan, initial conditions for each reactor, the inflow vector, and the number of reactors as input. The timespan is specified the same as in `ADM1sol`, the initial conditions are specified as a `Tuple` of vectors where each vector is the initial conditions for one of the reactors, and the inflow vector is specified the same as in `ADM1sol`. `MultiChamberSolution` returns the solutions to reactors as a `Tuple`, where the first element of the `Tuple` is the solution to reactor 1, the second is the solution to reactor 2, and so on.

The following code models three reactors in series. Each of the reactors has the same initial conditions.

```@repl
using ADM1jl

u0 = initialConditions(); # default initial conditions

IV = inflowvector_definition(); # default inflow vector

sols = MultiChamberSolution((0.0,200.0),(u0,u0,u0),IV,3); # solve the three reactors

sols[1] # solution to first reactor

sols[2] # solution to second reactor

sols[3] # solution to third reactor
```

### Variable Inflow

To specify variable inflow to the first reactor in the series, the input is much the same as in the `ADM1sol` case:

```@repl
using ADM1jl

u0 = initialConditions();

t = [i for i in 0.0:0.1:50.0];

IV_temp = inflowvector_definition();

IV = [IV_temp*(0.5*rand()+1.0) for i in 1:length(t)];

sols = MultiChamberSolution((0.0,50.0),(u0,u0,u0),IV,t,3);
```

## Plotting

```@repl
using ADM1jl

u0 = initialConditions(); # assigns the default initial conditions to u0

IV = inflowvector_definition(); # assigns the default inflow vector to IV

tspan = (0.0,200.0); # the solution will be computed from t=0.0 to t=200.0

sol, tSol = ADM1sol(tspan,u0,IV); # compute the solution

plotSols(sol)
```


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