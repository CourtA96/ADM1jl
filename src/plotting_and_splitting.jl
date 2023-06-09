function plotSols(sol)
   figure = plot(sol, vars=(0,1),reuse = false)
   for i = 2:size(sol)[1]
      plot!(sol,vars=(0,i))
   end
   display(figure)
end

function plotSols2(sol;titleText::String="Plots of Solutions")
   # I got the code to create the title from https://stackoverflow.com/questions/43066957/adding-global-title-to-plots-jl-subplots
   y = (ones(3))
   title = scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], text(titleText)),axis=false, grid=false, leg=false,size=(200,100),reuse=false)

   figure1temp = plot(sol,vars=1:20,layout=20,reuse=false,title = ["($i)" for j in 1:1, i in 1:20], titleloc = :right, titlefont = font(8),labels=:none)
   figure1 = plot(title,figure1temp,layout=grid(2,1,heights=[0.05,0.9]))
   display(figure1)

   sols = individualSolutions(sol)
   pltNums = [i for i in 21:size(sol)[1]]
   figure2temp = plot([(sol.t,sols[i]) for i in 21:35],layout=20,title = ["($i)" for j in 1:1, i in pltNums], titleloc = :right,reuse=false, titlefont = font(8),labels=:none)
   figure2 = plot(title,figure2temp,layout=grid(2,1,heights=[0.05,0.9]))
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
