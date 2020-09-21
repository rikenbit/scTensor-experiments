source("src/Functions.R")

# Parameter
# outfile <- "data/groundtruth/Human_Germline_Female.RData"
outfile <- commandArgs(trailingOnly=TRUE)[1]

# e.g. Human_Germline_Female
spl <- gsub("data/groundtruth/", "",
           gsub(".RData", "", outfile))

# Data loading
infile <- paste0("data/", spl, "/", spl, ".RData")
load(infile)
out <- scTensor:::.celltypemergedtensor(
    input=input,
    LR=LR,
    celltypes=celltypes,
    mergeas="mean",
    outerfunc="+")

# CCI1 : Female_Soma_3,Female_Soma_4 vs Female_FGC_2,Female_FGC_3,Female_FGC_4
# BMP2 (650) vs BMPR2 (659)
cci1 <- multicciInfo(
	lposition=c(6,4),
	rposition=c(5,7,8),
	ligandid="650",
	receptorid="659",
	out)

# CCI2 : Female_FGC_4 vs Female_Soma_1,Female_Soma_2,Female_Soma_3,Female_Soma_4
# JAG1 (182) vs NOTCH2 (4853)
cci2 <- multicciInfo(
	lposition=c(8),
	rposition=c(2,3,4,6),
	ligandid="182",
	receptorid="4853",
	out)

# True CaH and Color
trueCaH <- list(cci1$trueCaH, cci2$trueCaH)

# Number of celltypes e.g. 3
ncelltypes <- dim(out$tnsr)[1]

# CCI Information
cif <- list(
	nPair=dim(out$tnsr)[3],
	CCI1=list(
		LPattern=cci1$LPattern,
		RPattern=cci1$RPattern,
		nGene=1,
		fc="E2"),
	CCI2=list(
		LPattern=cci2$LPattern,
		RPattern=cci2$RPattern,
		nGene=1,
		fc="E2")
)

# Saving
save(trueCaH, ncelltypes, cif, spl, file=outfile)
