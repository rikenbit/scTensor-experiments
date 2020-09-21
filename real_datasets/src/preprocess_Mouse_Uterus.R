source("src/Functions.R")

########################################################
# 1. Gene expression matrix
########################################################

# Data loading
GSE118180 <- read.delim("data/Mouse_Uterus/GSM3320143_WT_Uterus_out_gene_exon_tagged.dge.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
row.names <- GSE118180[,1]
GSE118180 <- GSE118180[, 2:ncol(GSE118180)]
rownames(GSE118180) <- row.names

# Gene ID conversion
LefttoRight <- select(Mus.musculus,
  column=c("SYMBOL", "ENTREZID"),
  keytype="SYMBOL",
  keys=rownames(GSE118180))

GSE118180 <- convertToNCBIGeneID(
	GSE118180,
	rownames(GSE118180),
	LefttoRight)

# celltype label
label.GSE118180 <- read.delim("data/Mouse_Uterus/Uterus_metadata.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
label.GSE118180[,1] <- gsub("^.*\\.", "", label.GSE118180[, "Cell"])

common.barcode <- intersect(label.GSE118180[, "Cell"],
	colnames(GSE118180))
GSE118180 <- GSE118180[, common.barcode]
target <- unlist(sapply(common.barcode, function(x){
	which(x == label.GSE118180[, "Cell"])[1]
	}))
label.GSE118180 <- label.GSE118180[target, ]

# tSNE
res.tsne.GSE118180 <- label.GSE118180[, c("tSNE_1", "tSNE_2")]

# label
label.GSE118180 <- label.GSE118180[,4]
names(label.GSE118180) <- label.GSE118180
# Pink 1
label.GSE118180[
	which(label.GSE118180 == "Stromal")] <-
	brewer.pal(9, "PuRd")[6]
# Pink 2
label.GSE118180[
	which(label.GSE118180 == "StromalP")] <-
	brewer.pal(9, "PuRd")[8]
# Red
label.GSE118180[
	which(label.GSE118180 == "Endothelial")] <-
	brewer.pal(11, "Spectral")[1]
# Brown
label.GSE118180[
	which(label.GSE118180 == "Epithelial")] <-
	brewer.pal(11, "RdBu")[1]
# Green 1
label.GSE118180[
	which(label.GSE118180 == "Lymphocytes")] <-
	brewer.pal(11, "BrBG")[7]
# Green 2
label.GSE118180[
	which(label.GSE118180 == "Mesothelial")] <-
	brewer.pal(11, "BrBG")[8]
# Green 3
label.GSE118180[
	which(label.GSE118180 == "Lymphocytes Mesothelial")] <-
	brewer.pal(11, "BrBG")[9]
# Green 4
label.GSE118180[
	which(label.GSE118180 == "Myeloid")] <-
	brewer.pal(11, "BrBG")[10]
# Lymphatic Endothelial Cells : Green
label.GSE118180[
	which(label.GSE118180 == "LEC")] <-
	brewer.pal(11, "BrBG")[11]
# Blue
label.GSE118180[
	which(label.GSE118180 == "Myocytes")] <-
	brewer.pal(11, "RdBu")[10]
# Blue
label.GSE118180[
	which(label.GSE118180 == "Myocytes 2")] <-
	brewer.pal(11, "RdBu")[11]
# Purple
label.GSE118180[
	which(label.GSE118180 == "Pericytes")] <-
	brewer.pal(11, "Spectral")[11]
# Gray
label.GSE118180[
	which(label.GSE118180 == "11")] <-
	rgb(0.5,0.5,0.5)

# HVG
hvg.stats <- HVG(GSE118180)
target <- which(hvg.stats$pval < 1E-1)
GSE118180.HVG <- GSE118180[target, ]

geneid <- c("14560", "21809", "16001", "22339")
for(i in 1:length(geneid)){
	target <- which(rownames(GSE118180) == geneid[i])
	tmp <- GSE118180[target, ]
	GSE118180.HVG <- rbind(GSE118180.HVG, tmp)
}
nr1 <- nrow(GSE118180.HVG)
rownames(GSE118180.HVG)[(nr1-3):nr1] <- geneid

input <- GSE118180.HVG

########################################################
# 2. Celltype label vector
########################################################
celltypes <- label.GSE118180

########################################################
# 3. Ligand-Receptor pairs used in the original paper
########################################################
db <- c(
	"ENSEMBL_CELLPHONEDB", "ENSEMBL_SINGLECELLSIGNALR", "ENSEMBL_DLRP", "ENSEMBL_IUPHAR", "ENSEMBL_HPMR",
	"NCBI_CELLPHONEDB", "NCBI_SINGLECELLSIGNALR", "NCBI_DLRP", "NCBI_IUPHAR", "NCBI_HPMR")

tmp1 <- select(LRBase.Mmu.eg.db, columns=c("GENEID_L", "GENEID_R", "SOURCEDB"), keytype="GENEID_L", keys=rownames(input))
tmp2 <- select(LRBase.Mmu.eg.db, columns=c("GENEID_L", "GENEID_R", "SOURCEDB"), keytype="GENEID_R", keys=rownames(input))

tmp1 <- tmp1[which(tmp1$SOURCEDB %in% db), c("GENEID_L", "GENEID_R")]
tmp2 <- tmp2[which(tmp2$SOURCEDB %in% db), c("GENEID_L", "GENEID_R")]

LR <- unique(rbind(tmp1, tmp2))

########################################################
# Saving
########################################################
save(input, LR, celltypes,
	file="data/Mouse_Uterus/Mouse_Uterus.RData")
