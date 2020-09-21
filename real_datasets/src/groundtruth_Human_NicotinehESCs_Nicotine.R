source("src/Functions.R")

# Parameter
# outfile <- "data/groundtruth/Human_NicotinehESCs_Nicotine.RData"
# outfile <- commandArgs(trailingOnly=TRUE)[1]

# e.g. Human_NicotinehESCs_Nicotine
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

# CCI1 : Muscle vs Endothelial,Neural,USC
# GAS6 (2621) vs AXL (558)
cci1 <- multicciInfo(
	lposition=c(6),
	rposition=c(8,3,1),
	ligandid="2621",
	receptorid="558",
	out)

# CCI2 : Muscle vs Muscule
# COL4A1 (1282) vs CD47 (961)
cci2 <- multicciInfo(
	lposition=c(6),
	rposition=c(6),
	ligandid="1282",
	receptorid="961",
	out)

# CCI3 : Muscle vs USC
# EFNA1 (1942) vs EPHA1 (2041)
cci3 <- multicciInfo(
	lposition=c(6),
	rposition=c(1),
	ligandid="1942",
	receptorid="2041",
	out)

# CCI4 : Muscle vs Neural
# EFNA1 (1942) vs EPHA4 (2043)
cci4 <- multicciInfo(
	lposition=c(6),
	rposition=c(3),
	ligandid="1942",
	receptorid="2043",
	out)

# CCI5 : Endothelial,USC vs Muscle,USC
# FGF2 (2247) vs FGFR3 (2261)
cci5 <- multicciInfo(
	lposition=c(8,1),
	rposition=c(6,1),
	ligandid="2247",
	receptorid="2261",
	out)

# CCI6 : Epithelial,Liver,Neural vs Muscle,USC
# FGF8 (2253) vs FGFR3 (2261)
cci6 <- multicciInfo(
	lposition=c(2,4,3),
	rposition=c(6,1),
	ligandid="2253",
	receptorid="2261",
	out)

# CCI7 : Endothelial,USC vs Endothelial,Neural,USC
# COL6A1 (1291) vs ITGA6 (3655)
cci7 <- multicciInfo(
	lposition=c(8,1),
	rposition=c(8,3,1),
	ligandid="1291",
	receptorid="3655",
	out)

# CCI8 : Endothelial,USC vs Endothelial,Liver,Muscle,Neural,USC,Stromal
# COL1A1 (1277) vs ITGB1 (3688)
cci8 <- multicciInfo(
	lposition=c(8,1),
	rposition=c(8,4,6,3,1,5),
	ligandid="1277",
	receptorid="3688",
	out)

# CCI9 : Epithelial,Liver,Muscle,Neural,Stromal vs Endothelial,Epithelial,Liver,Muscle,Neural,USC,Stromal
# COL2A1 (1280) vs ITGB1 (3688)
cci9 <- multicciInfo(
	lposition=c(2,4,6,3,5),
	rposition=c(8,2,4,6,3,1,5),
	ligandid="1280",
	receptorid="3688",
	out)

# CCI10 : Muscle vs Endothelial,Epithelial,Liver,Muscle,Neural,USC,Stromal
# COL4A1 (1282) vs ITGB1 (3688)
cci10 <- multicciInfo(
	lposition=c(6),
	rposition=c(8,2,4,6,3,1,5),
	ligandid="1282",
	receptorid="3688",
	out)

# CCI11 : Endothelial,USC vs Endothelial,Epithelial,Liver,Muscle,Neural,USC,Stromal
# COL6A1 (1291) vs ITGB1 (3688)
cci11 <- multicciInfo(
	lposition=c(8,1),
	rposition=c(8,2,4,6,3,1,5),
	ligandid="1291",
	receptorid="3688",
	out)

# CCI12 : Endothelial,Epithelial,Liver,Muscle,Neural,USC,UDC,Stromal vs Endothelial,Epithelial,Liver,Muscle,Neural,USC,UDC,Stromal
# FBLN1 (2192) vs ITGB1 (3688)
cci12 <- multicciInfo(
	lposition=1:8,
	rposition=1:8,
	ligandid="2192",
	receptorid="3688",
	out)

# CCI13 : Endothelial vs Endothelial,Epithelial,Liver,Muscle,Neural,USC,Stromal
# PLAU (5328) vs ITGB1 (3688)
cci13 <- multicciInfo(
	lposition=c(8),
	rposition=c(8,2,4,6,3,1,5),
	ligandid="5328",
	receptorid="3688",
	out)

# CCI14 : USC vs Endothelial,Epithelial,Liver,Muscle,Neural,USC,Stromal
# CXCL12 (6387) vs ITGB1 (3688)
cci14 <- multicciInfo(
	lposition=c(1),
	rposition=c(8,2,4,6,3,1,5),
	ligandid="6387",
	receptorid="3688",
	out)

# CCI15 : Endothelial,USC vs Endothelial,Epithelial,Liver,Muscle,Neural,USC,Stromal
# SPP1 (6696) vs ITGB1 (3688)
cci15 <- multicciInfo(
	lposition=c(8,1),
	rposition=c(8,2,4,6,3,1,5),
	ligandid="6696",
	receptorid="3688",
	out)

# CCI16 : Epithelial,Liver,Neural vs Endothelial,Muscle,USC
# NRG1 (3084) vs ERBB3 (2065)
cci16 <- multicciInfo(
	lposition=c(2,4,3),
	rposition=c(8,6,1),
	ligandid="3084",
	receptorid="2065",
	out)

# CCI17 : Endothelial vs USC
# EFNA3 (1944) vs EPHA1 (2041)
cci17 <- multicciInfo(
	lposition=c(8),
	rposition=c(1),
	ligandid="1944",
	receptorid="2041",
	out)

# CCI18 : Endothelial,Epithelial,Liver,Muscle,Neural,USC,Stromal vs USC
# EFNA5 (1946) vs EPHA1 (2041)
cci18 <- multicciInfo(
	lposition=c(1,2,3,4,5,6,8),
	rposition=c(1),
	ligandid="1946",
	receptorid="2041",
	out)

# CCI19 : Epithelial,Liver,Muscle,Neural vs Neural
# EFNB1 (1947) vs EPHA4 (2043)
cci19 <- multicciInfo(
	lposition=c(2,4,6,3),
	rposition=c(3),
	ligandid="1947",
	receptorid="2043",
	out)

# CCI20 : Endothelial,Epithelial,Liver,Muscle,Neural,USC,Stromal vs Neural
# EFNA5 (1946) vs EPHA4 (2043)
cci20 <- multicciInfo(
	lposition=c(1,2,3,4,5,6,8),
	rposition=c(3),
	ligandid="1946",
	receptorid="2043",
	out)

# CCI21 : Endothelial vs Neural
# EFNA3 (1944) vs EPHA4 (2043)
cci21 <- multicciInfo(
	lposition=c(8),
	rposition=c(3),
	ligandid="1944",
	receptorid="2043",
	out)

# CCI22 : Epithelial,Liver,Neural vs Endothelial,Muscle,Neural,USC
# NRG1 (3084) vs ERBB2 (2064)
cci22 <- multicciInfo(
	lposition=c(2,4,3),
	rposition=c(8,6,3,1),
	ligandid="3084",
	receptorid="2064",
	out)

# True CaH and Color
trueCaH <- list(cci1$trueCaH, cci2$trueCaH, cci3$trueCaH, cci4$trueCaH, cci5$trueCaH, cci6$trueCaH, cci7$trueCaH, cci8$trueCaH, cci9$trueCaH, cci10$trueCaH, cci11$trueCaH, cci12$trueCaH, cci13$trueCaH, cci14$trueCaH, cci15$trueCaH, cci16$trueCaH, cci17$trueCaH, cci18$trueCaH, cci19$trueCaH, cci20$trueCaH, cci21$trueCaH, cci22$trueCaH)

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
		fc="E2"),
	CCI3=list(
		LPattern=cci3$LPattern,
		RPattern=cci3$RPattern,
		nGene=1,
		fc="E2"),
	CCI4=list(
		LPattern=cci4$LPattern,
		RPattern=cci4$RPattern,
		nGene=1,
		fc="E2"),
	CCI5=list(
		LPattern=cci5$LPattern,
		RPattern=cci5$RPattern,
		nGene=1,
		fc="E2"),
	CCI6=list(
		LPattern=cci6$LPattern,
		RPattern=cci6$RPattern,
		nGene=1,
		fc="E2"),
	CCI7=list(
		LPattern=cci7$LPattern,
		RPattern=cci7$RPattern,
		nGene=1,
		fc="E2"),
	CCI8=list(
		LPattern=cci8$LPattern,
		RPattern=cci8$RPattern,
		nGene=1,
		fc="E2"),
	CCI9=list(
		LPattern=cci9$LPattern,
		RPattern=cci9$RPattern,
		nGene=1,
		fc="E2"),
	CCI10=list(
		LPattern=cci10$LPattern,
		RPattern=cci10$RPattern,
		nGene=1,
		fc="E2"),
	CCI11=list(
		LPattern=cci11$LPattern,
		RPattern=cci11$RPattern,
		nGene=1,
		fc="E2"),
	CCI12=list(
		LPattern=cci12$LPattern,
		RPattern=cci12$RPattern,
		nGene=1,
		fc="E2"),
	CCI13=list(
		LPattern=cci13$LPattern,
		RPattern=cci13$RPattern,
		nGene=1,
		fc="E2"),
	CCI14=list(
		LPattern=cci14$LPattern,
		RPattern=cci14$RPattern,
		nGene=1,
		fc="E2"),
	CCI15=list(
		LPattern=cci15$LPattern,
		RPattern=cci15$RPattern,
		nGene=1,
		fc="E2"),
	CCI16=list(
		LPattern=cci16$LPattern,
		RPattern=cci16$RPattern,
		nGene=1,
		fc="E2"),
	CCI17=list(
		LPattern=cci17$LPattern,
		RPattern=cci17$RPattern,
		nGene=1,
		fc="E2"),
	CCI18=list(
		LPattern=cci18$LPattern,
		RPattern=cci18$RPattern,
		nGene=1,
		fc="E2"),
	CCI19=list(
		LPattern=cci19$LPattern,
		RPattern=cci19$RPattern,
		nGene=1,
		fc="E2"),
	CCI20=list(
		LPattern=cci20$LPattern,
		RPattern=cci20$RPattern,
		nGene=1,
		fc="E2"),
	CCI21=list(
		LPattern=cci21$LPattern,
		RPattern=cci21$RPattern,
		nGene=1,
		fc="E2"),
	CCI22=list(
		LPattern=cci22$LPattern,
		RPattern=cci22$RPattern,
		nGene=1,
		fc="E2")
)

# Saving
save(trueCaH, ncelltypes, cif, spl, file=outfile)
