source("src/Functions.R")

# Parameter
method <- commandArgs(trailingOnly=TRUE)[1]
E <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Plot
Plot_PR(method, E, outfile)
