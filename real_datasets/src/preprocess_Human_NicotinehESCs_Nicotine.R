source("src/Functions.R")

########################################################
# 1. Gene expression matrix
########################################################
# Data loading
# Control
GSE125416_Control <- readMM("data/Human_NicotinehESCs_Nicotine/GSM3573649_D_matrix.mtx")
GSE125416_Control <- as.matrix(GSE125416_Control)
# Column names
colnames(GSE125416_Control) <- unlist(read.delim("data/Human_NicotinehESCs_Nicotine/GSM3573649_D_barcodes.tsv", header=FALSE))
# Row names
rownames(GSE125416_Control) <- read.delim("data/Human_NicotinehESCs_Nicotine/GSM3573649_D_genes.tsv", header=FALSE)[,2]

# Nicotine
GSE125416_Nicotine <- readMM("data/Human_NicotinehESCs_Nicotine/GSM3573650_N_matrix.mtx")
GSE125416_Nicotine <- as.matrix(GSE125416_Nicotine)
# Column names
colnames(GSE125416_Nicotine) <- unlist(read.delim("data/Human_NicotinehESCs_Nicotine/GSM3573650_N_barcodes.tsv", header=FALSE))
# Row names
rownames(GSE125416_Nicotine) <- read.delim("data/Human_NicotinehESCs_Nicotine/GSM3573650_N_genes.tsv", header=FALSE)[,2]

# t-SNE and celltype label
res.tsne.GSE125416 <- read.delim("data/Human_NicotinehESCs_Nicotine/tSNE.txt",
	sep="\t", header=TRUE, stringsAsFactor=FALSE)

res.tsne.GSE125416_Control <- res.tsne.GSE125416[
	grep("^D_", res.tsne.GSE125416$X), c("tSNE_1", "tSNE_2")]
res.tsne.GSE125416_Nicotine <- res.tsne.GSE125416[
	grep("^N_", res.tsne.GSE125416$X), c("tSNE_1", "tSNE_2")]

# Sort
target_D <- paste0(gsub("^D_", "", res.tsne.GSE125416[grep("^D_", res.tsne.GSE125416$X), "X"]), "-1")
target_N <- paste0(gsub("^N_", "", res.tsne.GSE125416[grep("^N_", res.tsne.GSE125416$X), "X"]), "-1")
GSE125416_Control <- GSE125416_Control[, target_D]
GSE125416_Nicotine <- GSE125416_Nicotine[, target_N]

# GeneID conversion
LefttoRight <- select(Homo.sapiens,
  column=c("SYMBOL", "ENTREZID"),
  keytype="SYMBOL",
  keys=unique(c(rownames(GSE125416_Control),
  	rownames(GSE125416_Nicotine))))

GSE125416_Control <- convertRowID(GSE125416_Control,
	rownames(GSE125416_Control),
	LefttoRight)$output
GSE125416_Nicotine <- convertRowID(GSE125416_Nicotine,
	rownames(GSE125416_Nicotine),
	LefttoRight)$output

# HVG
hvg.stats <- HVG(cbind(GSE125416_Control, GSE125416_Nicotine))
target <- which(hvg.stats$pval < 1E-1)
GSE125416_Control.HVG <- GSE125416_Control[target, ]
GSE125416_Nicotine.HVG <- GSE125416_Nicotine[target, ]

geneid <- c("3146", "7099")
for(i in 1:length(geneid)){
	target1 <- which(rownames(GSE125416_Control) == geneid[i])
	target2 <- which(rownames(GSE125416_Nicotine) == geneid[i])
	tmp1 <- GSE125416_Control[target1, ]
	tmp2 <- GSE125416_Nicotine[target2, ]
	GSE125416_Control.HVG <- rbind(GSE125416_Control.HVG, tmp1)
	GSE125416_Nicotine.HVG <- rbind(GSE125416_Nicotine.HVG, tmp2)
}
nr1 <- nrow(GSE125416_Control.HVG)
nr2 <- nrow(GSE125416_Nicotine.HVG)
rownames(GSE125416_Control.HVG)[(nr1-1):nr1] <- geneid
rownames(GSE125416_Nicotine.HVG)[(nr2-1):nr2] <- geneid

input <- GSE125416_Nicotine.HVG

########################################################
# 2. Celltype label vector
########################################################
celltypes <- res.tsne.GSE125416[
	grep("^N_", res.tsne.GSE125416$X), "cellType"]
names(celltypes) <- celltypes
# Endothelial
celltypes[which(celltypes == "Endothelial")] <- ggdefault_cols(8)[1]
# Epithelial
celltypes[which(celltypes == "Epithelial")] <- ggdefault_cols(8)[2]
# Liver
celltypes[which(celltypes == "Liver")] <- ggdefault_cols(8)[3]
# Muscle
celltypes[which(celltypes == "Muscle")] <- ggdefault_cols(8)[4]
# Neural
celltypes[which(celltypes == "Neural")] <- ggdefault_cols(8)[5]
# USC
celltypes[which(celltypes == "USC")] <- ggdefault_cols(8)[6]
# UDC
celltypes[which(celltypes == "UDC")] <- ggdefault_cols(8)[7]
# Stromal
celltypes[which(celltypes == "Stromal")] <- ggdefault_cols(8)[8]

########################################################
# 3. Ligand-Receptor pairs used in the original paper
########################################################
# db <- c(
# 	"CELLPHONEDB", "SINGLECELLSIGNALR", "DLRP", "IUPHAR", "HPMR")
db <- c("DLRP", "IUPHAR", "HPMR")

# 2021/2/8書き換え部分
setAnnotationHubOption("CACHE", getwd())
ah <- AnnotationHub()
LRBase.Hsa.eg.db <- query(ah, c("LRBaseDb", "Homo sapiens", "v002"))[[1]]
LRBase.Hsa.eg.db <- LRBaseDb(LRBase.Hsa.eg.db)

tmp1 <- select(LRBase.Hsa.eg.db, columns=c("GENEID_L", "GENEID_R", "SOURCEDB"), keytype="GENEID_L", keys=rownames(input))
tmp2 <- select(LRBase.Hsa.eg.db, columns=c("GENEID_L", "GENEID_R", "SOURCEDB"), keytype="GENEID_R", keys=rownames(input))

tmp1 <- tmp1[which(tmp1$SOURCEDB %in% db), c("GENEID_L", "GENEID_R")]
tmp2 <- tmp2[which(tmp2$SOURCEDB %in% db), c("GENEID_L", "GENEID_R")]

LR <- unique(rbind(tmp1, tmp2))

########################################################
# Saving
########################################################
save(input, LR, celltypes,
		file="data/Human_NicotinehESCs_Nicotine/Human_NicotinehESCs_Nicotine.RData")