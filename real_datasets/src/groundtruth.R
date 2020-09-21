source("src/Functions.R")

# Parameter
# outfile <- "data/groundtruth/Human_Germline_Female.RData"
outfile <- commandArgs(trailingOnly=TRUE)[1]

# e.g. Mouse_Uterus
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

allgeneid <- c(
	unlist(lapply(out$pairname, function(x){
		strsplit(x, "_")[[1]][1]
	})),
	unlist(lapply(out$pairname, function(x){
		strsplit(x, "_")[[1]][2]
	})))
allgeneid <- unique(allgeneid)

# Setting
truepairs <- read.delim(paste0("data/truepairs/", spl, ".csv"), header=TRUE, stringsAsFactors=FALSE)
uniq.pairs <- unique(truepairs[,1:2])
ntypepairs <- nrow(uniq.pairs)

for(i in seq(ntypepairs)){
	print(i)
	print(uniq.pairs[i,1:2])
	if(length(grep("Mouse", spl)) == 1){
		ligandid <- as.character(select(Mus.musculus,
			columns="GENEID", keytype="SYMBOL",
			keys=uniq.pairs[i,1])[, "GENEID"])
		receptorid <- as.character(select(Mus.musculus,
			columns="GENEID", keytype="SYMBOL",
			keys=uniq.pairs[i,2])[, "GENEID"])
	}else{
		ligandid <- as.character(select(Homo.sapiens,
			columns="GENEID", keytype="SYMBOL",
			keys=uniq.pairs[i,1])[, "GENEID"])
		receptorid <- as.character(select(Homo.sapiens,
			columns="GENEID", keytype="SYMBOL",
			keys=uniq.pairs[i,2])[, "GENEID"])
	}
	ligandid <- intersect(ligandid, allgeneid)
	receptorid <- intersect(receptorid, allgeneid)
	target <- intersect(
		which(uniq.pairs[i,1] == truepairs[,1]),
		which(uniq.pairs[i,2] == truepairs[,2]))
	celltypesL <- unique(truepairs[target, 3])
	celltypesR <- unique(truepairs[target, 4])
	lposition <- sort(unique(unlist(lapply(celltypesL, function(x){
		which(x == unique(names(celltypes)))
	}))))
	rposition <- sort(unique(unlist(lapply(celltypesR, function(x){
		which(x == unique(names(celltypes)))
	}))))
	check1 <- length(lposition) == dim(out$tnsr)[1]
	check2 <- length(rposition) == dim(out$tnsr)[2]
	if(!check1 || !check2){
		# cci*
		eval(parse(text=paste(
			c("cci", i, " <- multicciInfo(lposition=lposition, rposition=rposition, ligandid=ligandid, receptorid=receptorid, out)"), collapse="")))
	}else{
		stop("All the celltypes are expressed!!!")
	}
	rm(ligandid)
	rm(receptorid)
}

# # True CaH
eval(parse(text=paste0(
	"trueCaH <- list(",
	paste0("cci", seq(ntypepairs), "$trueCaH", collapse=","),
	")")))

# Number of celltypes e.g. 3
ncelltypes <- dim(out$tnsr)[1]

# CCI Information
eval(parse(text=paste0("cif <- list(nPair=dim(out$tnsr)[3], ",
	paste0("CCI", seq(ntypepairs),
		" = list(LPattern=cci",
		seq(ntypepairs), "$LPattern,",
		" RPattern=cci",
		seq(ntypepairs), "$RPattern, nGene=1, fc='E2')"
		, collapse=","), ")")))

# Saving
save(trueCaH, ncelltypes, cif, spl, file=outfile)
