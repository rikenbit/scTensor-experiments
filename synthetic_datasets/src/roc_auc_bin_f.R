source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3] # ROC
outfile2 <- commandArgs(trailingOnly=TRUE)[4] # AUC
outfile3 <- commandArgs(trailingOnly=TRUE)[5] # BIN
outfile4 <- commandArgs(trailingOnly=TRUE)[6] # F
# 2022.01.25追加
outfile5 <- commandArgs(trailingOnly=TRUE)[7] # PRC
outfile6 <- commandArgs(trailingOnly=TRUE)[8] # AUCPR
outfile7 <- commandArgs(trailingOnly=TRUE)[9] # MCC
# 2022.02.24追加
outfile8 <- commandArgs(trailingOnly=TRUE)[10] # FPR
outfile9 <- commandArgs(trailingOnly=TRUE)[11] # FNR
# 2022.04.11追加
outfile10 <- commandArgs(trailingOnly=TRUE)[12] # PR

# Data loading
load(infile1)
load(infile2)

# Parameter
method <- strsplit(infile1, "/")[[1]][3]

# ROC / AUC
tensor.methods <- c(
	"labelpermutation_tensor", "labelpermutation2_tensor",
	"halpern_tensor", "cabelloaguilar_tensor")
pval.methods <- c(
	"labelpermutation", "labelpermutation2",
	"halpern", "cabelloaguilar")

if(method %in% tensor.methods){
    ROC_AUC_BIN_F(ncelltypes, cif, trueCaH, out,
    	outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9, outfile10,
    	pval=FALSE)
}
if(method %in% pval.methods){
    ROC_AUC_BIN_F(ncelltypes, cif, trueCaH, out,
    	outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9, outfile10,
    	pval=TRUE)
}
if(method == "previous_sctensor"){
    ROC_AUC_BIN_F_previous_sctensor(trueCaH, out,
    	outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9, outfile10)
}
if(method == "sctensor"){
    ROC_AUC_BIN_F_sctensor(trueCaH, out,
    	outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9, outfile10)
}
