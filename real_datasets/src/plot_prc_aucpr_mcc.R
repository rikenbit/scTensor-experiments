source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
infile4 <- commandArgs(trailingOnly=TRUE)[4]
outfile <- commandArgs(trailingOnly=TRUE)[5]

# Data loading
load(infile1)
load(infile2)
load(infile3)
load(infile4)

# Plot
Plot_PRC_AUCPR_MCC(trueCaH, prc, aucpr, mcc, outfile)
