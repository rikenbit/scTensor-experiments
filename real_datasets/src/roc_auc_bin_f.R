source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]
outfile3 <- commandArgs(trailingOnly=TRUE)[5]
outfile4 <- commandArgs(trailingOnly=TRUE)[6]
outfile5 <- commandArgs(trailingOnly=TRUE)[7]
outfile6 <- commandArgs(trailingOnly=TRUE)[8]
outfile7 <- commandArgs(trailingOnly=TRUE)[9]
outfile8 <- commandArgs(trailingOnly=TRUE)[10]
outfile9 <- commandArgs(trailingOnly=TRUE)[11]
outfile10 <- commandArgs(trailingOnly=TRUE)[12]

# Data loading
load(infile1)
load(infile2)

# Parameter
method <- strsplit(infile1, "/")[[1]][2]

# ROC / AUC
if(method %in% tensor.methods){
    ROC_AUC_BIN_F(ncelltypes, cif, trueCaH, out1,
    	outfile1, outfile2, outfile3, outfile4,
        outfile5, outfile6, outfile7, outfile8, outfile9, outfile10, pval=FALSE)
}
if(method %in% pval.methods){
    ROC_AUC_BIN_F(ncelltypes, cif, trueCaH, out2,
    	outfile1, outfile2, outfile3, outfile4,
        outfile5, outfile6, outfile7, outfile8, outfile9, outfile10, pval=TRUE)
}
if(method == "previous_sctensor"){
    ROC_AUC_BIN_F_previous_sctensor(trueCaH, out,
    	outfile1, outfile2, outfile3, outfile4,
        outfile5, outfile6, outfile7, outfile8, outfile9, outfile10)
}
if(method == "sctensor"){
    ROC_AUC_BIN_F_sctensor(trueCaH, out,
    	outfile1, outfile2, outfile3, outfile4,
        outfile5, outfile6, outfile7, outfile8, outfile9, outfile10)
}
