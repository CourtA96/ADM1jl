"""
# This file runs the julia  code for various inputs and outputs the result
# to a .xlsx file (ie. Microsoft Excel file)

# To read data back in:
# xf = XLSX.readxlsx("testADM1.xlsx")
# sheet = xf["Julia Data"]
# sheet[:] # returns all data contained the sheet as an array
# sheet[n+1,1] # returns the parameter number of the nth entry
# sheet[n+1,3] # returns the parameter value of the nth entry
# sheet[:][n+1,5:39] # returns the solution of the nth entry
"""

# If ADM1code hasn't already been included, do that
if !(:ADM1code in names(Main))
 include("ADM1code.jl")
end

# include modules for writing to excel files
import XLSX
using DataFrames
using DelimitedFiles
using DifferentialEquations

# initialize base u0 and inflow vector (IV) to work off of
u0standard = ADM1code.InitialConditions();
IVstandard = ADM1code.inflowvector_definition();

# Run the code once so it is loaded
ADM1code.MultiChamberSolutionExample((0.0,1.0),(u0standard,u0standard),IVstandard,2);

# Read data back in:
xf = XLSX.readxlsx("inputvalsRandInflow_15percent.xlsx")
inputDataSheet = xf[1]

global initialConditions = Array{Any}(undef,200)
global inflowVectors = Array{Any}(undef,200)

trialNums = [i for i in 1:400 if mod(i,2)==1]

for i in trialNums
    ind = Int((i+1)/2)
    IC_temp = inputDataSheet[:][i+1,3:37]
    IV_temp = inputDataSheet[:][i+2,3:37]
    global initialConditions[ind] = float.(IC_temp)
    global inflowVectors[ind] = float.(IV_temp)
end

XLSX.close(xf)

println("Starting")

# # For this to work, you need to have already created an excel spreadsheet
# # called "test.xlsx"
# XLSX.openxlsx("Julia Data Times 50.xlsx", mode="w") do xf # open the excel file
#     sheet = xf[1]
#     XLSX.rename!(sheet, "Random Inflow") # Rename the excel sheet to "Julia Data"
#     # Label the columns, the first colum gives the index of the parameter being
#     # changed, the third colum gives the value of the changed parameter,
#     # and the fifth column gives the final solution
#     temp = ["S_su", "S_aa","S_fa", "S_va", "S_bu", "S_pr", "S_ac", "S_h2",
#     "S_ch4", "S_IC", "S_IN", "S_I", "X_xc", "X_ch", "X_pr", "X_li", "X_su",
#     "X_aa", "X_fa", "X_c4", "X_pro", "X_ac", "X_h2", "X_I", "S_va_ion",
#     "S_bu_ion", "S_pro_ion", "S_ac_ion", "S_hco3_ion", "S_nh3", "S_cat",
#     "S_an", "S_gas_h2", "S_gas_ch4", "S_gas_co2", "time", "Initial Change", "Final Change", "Ratio of Final over Initial Change"]
#     sheet["A1"] = temp
# end

for i in 1:20 # Loop through the data

    # Reset u0 and IV each time a new variable is looked at
    u0 = initialConditions[i]
    IV = inflowVectors[i]

    # solve with the changed parameter
    sol,tS = ADM1code.MultiChamberSolutionExample((0.0,200.0),(u0,u0),IV,2);

    XLSX.openxlsx("RandInflow50with2Chambers.xlsx", mode="rw") do xf # open the excel file
        sheet = xf[1]
        # Save the information to the excel file

        sheet[string("A",i+1)] = sol[2][end]
        sheet[string("AJ",i+1)] = tS

        XLSX.close(xf)
    end

    println("Trial number ", i, " completed")
end
