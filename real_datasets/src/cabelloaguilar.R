source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1] 
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Data loading
load(infile)

# Label Permutation Tensor
out1 <- scTensor:::.celltypemergedtensor.ca(
    input=input,
    LR=LR,
    celltypes=celltypes,
    mergeas="mean",
    outerfunc="*")

# Label Permutation
out2 <- scTensor:::.cellCellDecomp.CabelloAguilar(
    input=input,
    LR=LR,
    celltypes=celltypes,
    ranks=c(3,3,3),
    rank=3,
    centering=TRUE,
    mergeas="mean",
    outerfunc="*",
    comb="random",
    num.sampling=100,
    num.perm=1000,
    decomp=TRUE,
    thr1=log2(5),
    thr2=25)

# Save
save(out1, file=outfile1)
save(out2, file=outfile2)
