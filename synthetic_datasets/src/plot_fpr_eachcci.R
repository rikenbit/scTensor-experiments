source("src/Functions.R")

# Parameter
method <- commandArgs(trailingOnly=TRUE)[1]
E <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Plot
Plot_FPR_eachCCI(method, E, outfile)
