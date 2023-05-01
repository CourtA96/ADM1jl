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

# If ADM1codeDAE hasn't already been included, do that
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
ADM1code.ExampleSol((0.0,1.0),u0standard,IVstandard);

# Vector that hold the factors by which each each parameter will be multiplied
changes = [0.01,0.1,10,100]

println("Starting")

# For this to work, you need to have already created an excel spreadsheet
# called "test.xlsx"
XLSX.openxlsx("test.xlsx", mode="w") do xf # open the excel file
    sheet = xf[1]
    XLSX.rename!(sheet, "Julia Data") # Rename the excel sheet to "Julia Data"
    # Label the columns, the first colum gives the index of the parameter being
    # changed, the third colum gives the value of the changed parameter,
    # and the fifth column gives the final solution
    sheet["A1"] = "Param Number"
    sheet["C1"] = "Param Value"
    sheet["E1"] = "Final Solution"

    global cell = 2
    XLSX.close(xf)
end

for i in 1:35 # Loop through the variables

    # Reset u0 and IV each time a new variable is looked at
    u0 = [j for j in u0standard]
    IV = [j for j in IVstandard]

    for change in changes # loop over the multiplication factors in changes
        print("Parameter ",  i, " x ",change,": ")

        param = string(i)

        # multiply the value of parameter i by the factor
        u0[i] = u0standard[i]*change
        IV[i] = IVstandard[i]*change

        # solve with the changed parameter
        sol = ADM1code.ExampleSol((0.0,200.0),u0,IV);

        XLSX.openxlsx("test.xlsx", mode="rw") do xf # open the excel file
            sheet = xf[1]
            # Save the information to the excel file
            sheet[string("A",cell)] = i
            sheet[string("C",cell)] = u0[i]
            sheet[string("E",cell)] = sol[end]

            XLSX.close(xf)
        end

        print("completed \n")

        global cell += 1
    end
    println("*** Parameter ", i, " completed ***")
end

println("Done")
