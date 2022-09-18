source("src/Functions.R")

# Parameter
method <- commandArgs(trailingOnly=TRUE)[1]
sample <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Plot
Plot_LR(method, sample, outfile)
