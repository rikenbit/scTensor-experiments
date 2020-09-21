source("src/Functions.R")

# Parameter
# infile1 = "output/sctensor/Human_NicotinehESCs_Nicotine.RData"
# infile2 = "data/groundtruth/Human_NicotinehESCs_Nicotine.RData"
# outfile1 = "output/sctensor/ROC/Human_NicotinehESCs_Nicotine.RData"
# outfile2 = "output/sctensor/AUC/Human_NicotinehESCs_Nicotine.RData"
# outfile3 = "output/sctensor/BIN/Human_NicotinehESCs_Nicotine.RData"
# outfile4 = "output/sctensor/F/Human_NicotinehESCs_Nicotine.RData"

infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]
outfile3 <- commandArgs(trailingOnly=TRUE)[5]
outfile4 <- commandArgs(trailingOnly=TRUE)[6]

# Data loading
load(infile1)
load(infile2)

# Parameter
method <- strsplit(infile1, "/")[[1]][2]

# ROC / AUC
if(method %in% tensor.methods){
    ROC_AUC_BIN_F(ncelltypes, cif, trueCaH, out1,
    	outfile1, outfile2, outfile3, outfile4, pval=FALSE)
}
if(method %in% pval.methods){
    ROC_AUC_BIN_F(ncelltypes, cif, trueCaH, out2,
    	outfile1, outfile2, outfile3, outfile4, pval=TRUE)
}
if(method == "previous_sctensor"){
    ROC_AUC_BIN_F_previous_sctensor(trueCaH, out,
    	outfile1, outfile2, outfile3, outfile4)
}
if(method == "sctensor"){
    ROC_AUC_BIN_F_sctensor(trueCaH, out,
    	outfile1, outfile2, outfile3, outfile4)
}
