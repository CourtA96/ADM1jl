# List of Functions
```@meta
CurrentModule = ADM1jl
```

## Single Reactor Solvers

```@docs
ADM1sol(tspan::Tuple,u0::Vector,IV::Vector{Float64};alg = Rodas4P(), tols=1e-4,tMax = 300.0)
```

```@docs
ADM1sol(tspan::Tuple,u0::Vector,IV::Vector{Vector{Float64}},IVtimes::Vector{Float64};alg = Rodas4P(), tols=1e-4,tMax = 300.0)
```

## Multiple Reactor Solvers

```@docs
MultiChamberSolution(tspan::Tuple,u0::Tuple,IV::Vector,nChambers::Int64;alg = Rodas4P(), tols=1e-4,tMax = 300.0)
```




```@docs
RHSfun(du,u,p,t)
```

```@docs
reactionrates(bp,rp,php,pressures,sx,NREAC::Int)
```

```@docs
pressureOfGasses(sx,php,rp)
```

```@docs
monod(u, k)
```