source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]
spl <- gsub("data/groundtruth/", "",
           gsub(".RData", "", outfile))

# Ground Truth CaH and Color
tmp <- groundTruth(spl)
trueCaH <- tmp$trueCaH
ncelltypes <- tmp$ncelltypes
cif <- tmp$cif
spl <- tmp$spl

# Saving
save(trueCaH, ncelltypes, cif, spl, file=outfile)
