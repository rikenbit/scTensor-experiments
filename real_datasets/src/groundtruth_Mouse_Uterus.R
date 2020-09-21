source("src/Functions.R")

# Parameter
# outfile <- "data/groundtruth/Mouse_Uterus.RData"
outfile <- commandArgs(trailingOnly=TRUE)[1]

# e.g. Mouse_Uterus
spl <- gsub("data/groundtruth/", "",
           gsub(".RData", "", outfile))

# True pairs
truepairs <- read.delim("data/truepairs/Mouse_Uterus.csv")

# Data loading
infile <- paste0("data/", spl, "/", spl, ".RData")
load(infile)
out <- scTensor:::.celltypemergedtensor(
    input=input,
    LR=LR,
    celltypes=celltypes,
    mergeas="mean",
    outerfunc="+")

# CCI1 : Stromal vs Epithelial
# Rspo3 (72780) vs Lgr5 (14160)
cci1 <- multicciInfo(
	lposition=c(1),
	rposition=c(2),
	ligandid="72780",
	receptorid="14160",
	out)

# CCI2 : Stromal vs Mesothelial,Myocytes,Myocytes 2,Epithelial
# Igf1 (16000) vs Igf1r (16001)
cci2 <- multicciInfo(
	lposition=c(1),
	rposition=c(11,4,5,2),
	ligandid="16000",
	receptorid="16001",
	out)

# True CaH
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
