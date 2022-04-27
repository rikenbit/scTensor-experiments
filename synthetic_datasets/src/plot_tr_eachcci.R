source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]

# Plot
Plot_TR_eachCCI(outfile)
