function plotSols(sol;titleText::String="Plots of Solutions")
   # I got the code to create the title from https://stackoverflow.com/questions/43066957/adding-global-title-to-plots-jl-subplots
   y = (ones(3))
   title1 = scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], text(string(titleText, " (1 of 2)"),20,"Helvetica")),axis=false, grid=false, leg=false,size=(200,100), labels=:none,reuse=false)
   title2 = scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], text(string(titleText, " (2 of 2)"),20,"Helvetica")),axis=false, grid=false, leg=false,size=(200,200), labels=:none,reuse=false)

   sols = individualSolutions(sol)

   figure1temp = plot([(sol.t,sols[i]) for i in 1:20],layout=20,reuse=false,title = ["($i)" for j in 1:1, i in 1:20], titleloc = :left, titlefont = font(8,"Helvetica Bold"), tickfontsize = 8, labels=:none, size = (1000, 700))
   figure1 = plot(title1,figure1temp,layout=grid(2,1,heights=[0.05,0.9]))
   display(figure1)

   pltNums = [i for i in 21:size(sol)[1]]
   figure2temp = plot([(sol.t,sols[i]) for i in 21:35],layout=(3,5),title = ["($i)" for j in 1:1, i in pltNums], titleloc = :left,reuse=false, titlefont = font(8,"Helvetica Bold"), tickfontsize = 8, labels=:none, size = (1000, 550))
   figure2 = plot(title2,figure2temp,layout=grid(2,1,heights=[0.05,0.9]))
   display(figure2)
end

function individualSolutions(sol)
   """
   Returns a 1D array containing vectors of each solution over time.

      i.e.
      A[i] is a vector containing the values of the ith solution for all
      times.

      So, if we say t is a vector containing the timesteps, then
         plot(t,A[i])
      will plot the ith solution vs time.

   """

    A = [[sol[i][j] for i in 1:length(sol.t)] for j in 1:size(sol)[1]]
    return A

end
