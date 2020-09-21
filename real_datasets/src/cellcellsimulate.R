source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]
spl <- gsub("data/", "", gsub(".RData", "", outfile))
cif <- cciInfo[[spl]]
ncl <- nCells[[spl]]

# Setting
params <- newCCSParams()
setParam(params, "nGene") <- 10000
setParam(params, "nCell") <- ncl
setParam(params, "cciInfo") <- cif

# Generating simulation data
out <- try(cellCellSimulate(params))
input <- out$input
LR <- out$LR
celltypes <- out$celltypes
LR_CCI <-out$LR_CCI

# Saving
save(input, LR, celltypes, LR_CCI, file=outfile)
