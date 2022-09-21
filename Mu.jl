include("Control.jl")

using Dates

df = DateFormat("dd-u-yyyy HH:MM:SS.sss")
println("Program started at " * Dates.format(Dates.now(), df))

Control.calculate()

println("Program ended at " * Dates.format(Dates.now(), df))
