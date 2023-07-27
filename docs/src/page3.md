# List of Functions
```@meta
CurrentModule = ADM1jl
```

## Single Reactor Solvers

```@docs
ADM1sol
```

## Multiple Reactor Solvers

```@docs
MultiChamberSolution
```
## Other Functions

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