source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3] # ROC
outfile2 <- commandArgs(trailingOnly=TRUE)[4] # AUC
outfile3 <- commandArgs(trailingOnly=TRUE)[5] # BIN
outfile4 <- commandArgs(trailingOnly=TRUE)[6] # F

# infile1 = "output/E10/halpern/30Celltypes_5CCIPatterns_ManytoMany.RData"
# infile2 = "data/groundtruth/30Celltypes_5CCIPatterns_ManytoMany.RData"

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
    	outfile1, outfile2, outfile3, outfile4, pval=FALSE)
}
if(method %in% pval.methods){
    ROC_AUC_BIN_F(ncelltypes, cif, trueCaH, out,
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
