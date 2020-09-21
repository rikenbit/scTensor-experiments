source("src/Functions.R")

# Parameter
E <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

spl <- gsub(paste0("data/", E, "/"),
	"",
	gsub(".RData", "", outfile))
cif <- cciInfo[[E]][[spl]]
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

# Saving
save(input, LR, celltypes, file=outfile)
