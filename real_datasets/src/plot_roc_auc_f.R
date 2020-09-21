source("src/Functions.R")

# Parameter
# infile1 = "output/sctensor/ROC/Human_NicotinehESCs_Nicotine.RData"
# infile2 = "output/sctensor/AUC/Human_NicotinehESCs_Nicotine.RData"
# infile3 = "output/sctensor/F/Human_NicotinehESCs_Nicotine.RData"
# outfile = "plot/ROC/sctensor/Human_NicotinehESCs_Nicotine.png"

infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
outfile <- commandArgs(trailingOnly=TRUE)[4]

# Data loading
load(infile1)
load(infile2)
load(infile3)

# Plot
Plot_ROC_AUC_F(roc, auc, f, outfile)
