source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Data Loading
load(infile)

# Normalization
input <- log10(scTensor:::.CPMED(input) + 1)

# Rank estimation（cellCellRanksの中身）
tnsr <- scTensor:::.cellCellDecomp.Third(
    input,
    LR,
    celltypes,
    ranks=c(1,1,1),
    rank=1,
    centering=TRUE,
    mergeas="mean",
    outerfunc="*",
    comb="random",
    num.sampling=100,
    num.perm=1000,
    decomp=FALSE,
    thr1=log2(5),
    thr2=25,
    thr3=0.95,
    L1_A=0,
    L2_A=0,
    verbose=TRUE)$cellcelllrpairpattern

l <- length(unique(names(celltypes)))

l1 <- min(dim(tnsr)[1], dim(tnsr)[2]*dim(tnsr)[3])
l2 <- min(dim(tnsr)[2], dim(tnsr)[3]*dim(tnsr)[1])
l3 <- min(10, dim(tnsr)[3], dim(tnsr)[1]*dim(tnsr)[2])

out1 <- NMF(cs_unfold(tnsr, m=1)@data, runtime=5, rank.method="rss", J=1:l1)
out2 <- NMF(cs_unfold(tnsr, m=2)@data, runtime=5, rank.method="rss", J=1:l2)
out3 <- NMF(cs_unfold(tnsr, m=3)@data, runtime=5, rank.method="rss", J=1:l3)

rss1 <- unlist(lapply(seq(l1), function(x, out1){
    eval(parse(text=paste0("mean(out1$Trial$Rank", x, "$original)")))
}, out1=out1))
rss2 <- unlist(lapply(seq(l2), function(x, out2){
    eval(parse(text=paste0("mean(out2$Trial$Rank", x, "$original)")))
}, out2=out2))
rss3 <- unlist(lapply(seq(l3), function(x, out3){
    eval(parse(text=paste0("mean(out3$Trial$Rank", x, "$original)")))
}, out3=out3))

rank1 <- min(which((max(rss1) - rss1) / (max(rss1) - min(rss1)) > 0.8))
rank2 <- min(which((max(rss2) - rss2) / (max(rss2) - min(rss2)) > 0.8))
rank3 <- min(which((max(rss3) - rss3) / (max(rss3) - min(rss3)) > 0.8))

selected = c(rank1, rank2, rank3)
print(selected)

# Perform scTensor
out <- scTensor:::.cellCellDecomp.Third(
    input=input,
    LR=LR,
    celltypes=celltypes,
    ranks=selected,
    rank=3,
    centering=TRUE,
    mergeas="mean",
    outerfunc="*",
    comb="random",
    num.sampling=100,
    num.perm=1000,
    decomp=TRUE,
    thr1=log2(5),
    thr2=25,
    thr3=0.95,
    L1_A=0,
    L2_A=0,
    verbose=TRUE)

# Save
save(out, file=outfile)