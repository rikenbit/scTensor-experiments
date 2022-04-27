source("src/Functions.R")

# Parameter
E <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Plot
Plot_PR_Merge(E, outfile)
