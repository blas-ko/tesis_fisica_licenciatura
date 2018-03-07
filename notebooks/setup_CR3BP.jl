using TaylorSeries, TaylorIntegration, Distributions
using Plots, LaTeXStrings
pyplot()

println("Using TaylorSeries, TaylorIntegration, Distributions, DataFrames, Plots, LaTeXStrings packages... \n")

## Custom libraries ###
include("JTFunctions.jl")
using .JTFunctions
include("ThreeBody.jl")
using .ThreeBody
import ThreeBody: μ

include("params_CR3BP.jl")
