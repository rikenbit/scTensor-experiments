source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1] 
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Data loading
load(infile)

# Label Permutation
out <- scTensor:::.cellCellDecomp.LabelPerm.LR(
	input=input,
	LR=LR,
	celltypes=celltypes,
	ranks=c(3,3,3),
    rank=3,
    centering=TRUE,
    mergeas="mean",
    outerfunc="+",
    comb="random",
    num.sampling=100,
    num.perm=1000,
    decomp=TRUE,
    thr1=log2(5),
    thr2=25)

# Save
save(out, file=outfile)
