export plotSols
"""
    plotSols(sol;<keyword arguments>)

Plot the solutions returned by ADM1sol. The plots are split between two figures and are 
displayed by default. However, due to an error in Julia's Plots package, the plot's windows 
may overwrite eachother in some circumstances. If this occurs, use the keyword arguments 
to view both plots manually.

# Arguments
- `sol::Vector`: `ODESolution` returned by ADM1 sol

# Optional Arguments
- `titleText = "Plot of Solutions"`: The title of the plots
- `displayPlots = true`: Boolean to display the plots or not.
- `savePNG = false`: Boolean to save pngs of the plots to the working directory as titleText(1 of 2).png and titleText(2 of 2).png
- `returnPlots = false`: Boolean to return the plots as a tuple so they can be displayed manually

# Examples
```julia-repl
julia> u0 = initialConditions();

julia> IV = inflowvector_definition();

julia> sol, tSol = ADM1sol((0.0,200.0),u0,IV); # compute the solution

julia> plotSols(sol); # display the plots
```

```julia-repl
julia> u0 = initialConditions();

julia> IV = inflowvector_definition();

julia> sol, tSol = ADM1sol((0.0,200.0),u0,IV); # compute the solution

julia> plt1,plt2 = plotSols(sol,displayPlots=false,returnPlots=true); # return the plots

julia> display(plt1); # display the first plot

julia> display(plt2); # display the second plot
```
"""
function plotSols(sol;titleText::String="Plots of Solutions",displayPlots=true,savePNG=false,returnPlots=false)
   # I got the code to create the title from https://stackoverflow.com/questions/43066957/adding-global-title-to-plots-jl-subplots
    y = (ones(3))
    title1 = scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], text(string(titleText, " (1 of 2)"),20,"Helvetica")),axis=false, grid=false, leg=false,size=(200,100), labels=:none,reuse=false)
    title2 = scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], text(string(titleText, " (2 of 2)"),20,"Helvetica")),axis=false, grid=false, leg=false,size=(200,200), labels=:none,reuse=false)

    sols = individualSolutions(sol)

    figure1temp = plot([(sol.t,sols[i]) for i in 1:20],layout=20,reuse=false,title = ["($i)" for j in 1:1, i in 1:20], titleloc = :left, titlefont = font(8,"Helvetica Bold"), tickfontsize = 8, labels=:none, size = (1000, 700))
    figure1 = plot(title1,figure1temp,layout=grid(2,1,heights=[0.05,0.9]),reuse=false)

    pltNums = [i for i in 21:size(sol)[1]]
    figure2temp = plot([(sol.t,sols[i]) for i in 21:35],layout=(3,5),title = ["($i)" for j in 1:1, i in pltNums], titleloc = :left,reuse=false, titlefont = font(8,"Helvetica Bold"), tickfontsize = 8, labels=:none, size = (1000, 550))
    figure2 = plot(title2,figure2temp,layout=grid(2,1,heights=[0.05,0.9]),reuse=false)

    if savePNG == true
        png(figure1,string(titleText," (1 of 2)"))
        png(figure2,string(titleText," (2 of 2)"))
    end

    if displayPlots == true
        display(figure1)
        display(figure2)
    end

    if returnPlots == true
        return figure1,figure2
    end
end

export individualSolutions
"""
    individualSolutions(sol)

Rearranges the `ODESolution` returned by `ADM1sol` into  a 1D array containing 
vectors of each solution over time.

i.e.
A[i] is a vector containing the values of the ith solution for all times.

So, if we say t is a vector containing the timesteps, then
   `plot(t,A[i])`
will plot the ith solution vs time.

# Arguments
- `sol::Vector`: `ODESolution` returned by ADM1 sol

# Examples
```julia-repl
julia> u0 = initialConditions();

julia> IV = inflowvector_definition();

julia> sol, tSol = ADM1sol((0.0,200.0),u0,IV); # compute the solution

julia> indSols = individualSolutions(sol)
35-element Vector{Vector{Float64}}:
 [0.012, 0.011999986305998722, 0.011999977102201841, 0.01199997267285494, 0.011999969386113186, 0.011999967149464337, 0.011999965091000644, 0.011999963258166062, 0.01199996178749697, 0.011999961662742068  …  0.01195750330094027, 0.011957024281388772, 0.01195654102848209, 0.011955980109704207, 0.011955371939386264, 0.011955036272014314, 0.011954879037458292, 0.011954836213455822, 0.011954830158706282, 0.011954829810563401]
 [0.0053, 0.005299528080754454, 0.005299211849244109, 0.005299059933102669, 0.005298947319276687, 0.005298870740171685, 0.005298800301401363, 0.005298737615364131, 0.005298687337685513, 0.005298683073531872  …  0.00531290262808903, 0.005313145376836529, 0.005313458937485294, 0.00531386295705003, 0.0053143158605383585, 0.0053145727297110815, 0.005314697865111298, 0.00531473416356305, 0.0053147398286552405, 0.0053147401832121036]
 ⋮
 [1.63, 1.6300010043803432, 1.630003893522008, 1.6300064891566448, 1.6300090584027846, 1.630011157417349, 1.630013359837308, 1.6300155516628245, 1.6300174749013072, 1.6300176449752604  …  1.6258740229661408, 1.6256583127128807, 1.6256250660754588, 1.625621749084325, 1.625620443499663, 1.6256197935003422, 1.6256194976733636, 1.6256194164919724, 1.6256194045086088, 1.6256194037384155]
 [0.014, 0.013987846051651676, 0.01395723325778608, 0.013936637777305555, 0.013919032066060481, 0.013905958298756526, 0.01389316476641871, 0.013881173393149072, 0.013871152739877324, 0.013870287064082722  …  0.01414633009687048, 0.014149764848711631, 0.01415028486767014, 0.014150328062007254, 0.01415033918241631, 0.014150344102072282, 0.014150346161797715, 0.014150346675181056, 0.014150346741486993, 0.014150346745679485]

 julia> indSols[1] # Solution vector for the first state variable (S_su)
 115-element Vector{Float64}:
 0.012
 0.011999986305998722
 ⋮
 0.011954830158706282
 0.011954829810563401
```
"""
function individualSolutions(sol)

    A = [[sol[i][j] for i in 1:length(sol.t)] for j in 1:size(sol)[1]]
    return A

end
