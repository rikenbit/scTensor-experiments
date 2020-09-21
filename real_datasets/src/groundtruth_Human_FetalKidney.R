source("src/Functions.R")

# Parameter
# outfile <- "data/groundtruth/Human_FetalKidney.RData"
outfile <- commandArgs(trailingOnly=TRUE)[1]

# e.g. Human_FetalKidney
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


# CCI1 : Podocyte, Collecting duct, Endothelial vs Cap mesenchyme, Nephron progenitor, Loop of Henle / Distal, Connecting tubule / Distal
# EFNB2 (1948) vs EPHB3 (2049)
cci1 <- multicciInfo(
	lposition=c(10,7,4),
	rposition=c(1,6,5,11),
	ligandid="1948",
	receptorid="2049",
	 out)

# CCI2 : Proximal tubule vs Endothelial
# JAG1 (182) vs NOTCH4 (4855)
cci2 <- multicciInfo(
	lposition=8,
	rposition=4,
	ligandid="1948",
	receptorid="2049",
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
