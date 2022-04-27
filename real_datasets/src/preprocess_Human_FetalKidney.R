source("src/Functions.R")

########################################################
# 1. Gene expression matrix
########################################################

# Data Loading
tmp <- read.delim("data/Human_FetalKidney/GSM3300889_FetalKidney_Run1865_87days.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactor=FALSE)
tmp2 <- read.delim("data/Human_FetalKidney/GSM2935201_FetalKidney_Run1751_105days.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactor=FALSE)
tmp3 <- read.delim("data/Human_FetalKidney/GSM2935200_FetalKidney_Run1785_110days.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactor=FALSE)
tmp4 <- read.delim("data/Human_FetalKidney/GSM2935199_FetalKidney_Run1824_115days.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactor=FALSE)
tmp5 <- read.delim("data/Human_FetalKidney/GSM3300888_FetalKidney_Run1865_132days.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactor=FALSE)

row.names <- intersect(rownames(tmp),
				intersect(rownames(tmp2),
					intersect(rownames(tmp3),
						intersect(rownames(tmp4),
							rownames(tmp5)))))

GSE109205 <- cbind(tmp[row.names,],
				tmp2[row.names,],
					tmp3[row.names,],
						tmp4[row.names,],
							tmp5[row.names,])

# t-SNE
blacklist <- c("ACGACATGATAG", "CCCCAAACGGTC", "CCCCCTAGCCGA")
res.tsne.GSE109205 <- read.delim("data/Human_FetalKidney/FK_TSNE_co-ordinates.dms", sep="\t", header=TRUE, stringsAsFactor=FALSE)
res.tsne.GSE109205 <- res.tsne.GSE109205[
	setdiff(1:nrow(res.tsne.GSE109205),
		as.vector(sapply(blacklist,
			function(x){grep(x, res.tsne.GSE109205[,"X"])}))), ]
tmp <- res.tsne.GSE109205[,"X"]
res.tsne.GSE109205 <- res.tsne.GSE109205[, c("tSNE_1", "tSNE_2")]
rownames(res.tsne.GSE109205) <- gsub("^r.*__", "", tmp)

# celltype label
clusternames <- read.delim("data/Human_FetalKidney/FK_FirstClusterNames.txt", stringsAsFactor=FALSE)
label.GSE109205 <- read.delim("data/Human_FetalKidney/FK_Barcodes_FirstClustering.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
target <- setdiff(1:nrow(label.GSE109205),
	unlist(sapply(11:18, function(x){which(x == label.GSE109205)})))
barcode <- rownames(label.GSE109205)[target]
barcode <- gsub("^r.*__", "", barcode)
label.GSE109205 <- label.GSE109205[target, ]
names(label.GSE109205) <- barcode

# Sort
target <- intersect(unique(names(label.GSE109205)),
	unique(rownames(res.tsne.GSE109205)))
label.GSE109205 <- label.GSE109205[target]
GSE109205 <- GSE109205[, target]

# Gene ID conversion
LefttoRight <- select(Homo.sapiens,
  column=c("SYMBOL", "ENTREZID"),
  keytype="SYMBOL",
  keys=rownames(GSE109205))

GSE109205 <- convertRowID(
	GSE109205,
	rownames(GSE109205),
	LefttoRight)$output

# Color
for(i in 0:10){
	label.GSE109205[which(label.GSE109205 == i)] <-
		clusternames[which(clusternames$Cluster == i), "Celltype"]
}
names(label.GSE109205) <- label.GSE109205

# 0: Nephron progenitor
label.GSE109205[which(label.GSE109205 == "Nephron progenitor")] <- ggdefault_cols(11)[1]
# 1: Cap mesenchyme
label.GSE109205[which(label.GSE109205 == "Cap mesenchyme")] <- ggdefault_cols(11)[2]
# 2: Stroma
label.GSE109205[which(label.GSE109205 == "Stroma")] <- ggdefault_cols(11)[3]
# 3: Loop of Henle / Distal
label.GSE109205[which(label.GSE109205 == "Loop of Henle / Distal")] <- ggdefault_cols(11)[4]
# 4: Collecting duct
label.GSE109205[which(label.GSE109205 == "Collecting duct")] <- ggdefault_cols(11)[5]
# 5: Endothelial
label.GSE109205[which(label.GSE109205 == "Endothelial")] <- ggdefault_cols(11)[6]
# 6: Proximal tubule
label.GSE109205[which(label.GSE109205 == "Proximal tubule")] <- ggdefault_cols(11)[7]
# 7: Podocyte
label.GSE109205[which(label.GSE109205 == "Podocyte")] <- ggdefault_cols(11)[8]
# 8: Connecting tubule / Distal
label.GSE109205[which(label.GSE109205 == "Connecting tubule / Distal")] <- ggdefault_cols(11)[9]
# 9: Unspecified/RBC
label.GSE109205[which(label.GSE109205 == "Unspecified/RBC")] <- ggdefault_cols(11)[10]
# 10: Immune
label.GSE109205[which(label.GSE109205 == "Immune")] <- ggdefault_cols(11)[11]

# HVG
hvg.stats <- HVG(GSE109205)
target <- which(hvg.stats$pval < 1E-1)
GSE109205.HVG <- GSE109205[target, ]

geneid <- c("7048", "4855")
for(i in 1:length(geneid)){
	target <- which(rownames(GSE109205) == geneid[i])
	tmp <- GSE109205[target, ]
	GSE109205.HVG <- rbind(GSE109205.HVG, tmp)
}

input <- GSE109205.HVG

########################################################
# 2. Celltype label vector
########################################################
celltypes <- label.GSE109205

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

tmp1 <- select(LRBase.Hsa.eg.db,
	columns=c("GENEID_L", "GENEID_R", "SOURCEDB"),
	keytype="GENEID_L", keys=rownames(input))
tmp2 <- select(LRBase.Hsa.eg.db,
	columns=c("GENEID_L", "GENEID_R", "SOURCEDB"),
	keytype="GENEID_R", keys=rownames(input))

tmp1 <- tmp1[which(tmp1$SOURCEDB %in% db), c("GENEID_L", "GENEID_R")]
tmp2 <- tmp2[which(tmp2$SOURCEDB %in% db), c("GENEID_L", "GENEID_R")]

LR <- unique(rbind(tmp1, tmp2))

########################################################
# Saving
########################################################
save(input, LR, celltypes,
	file="data/Human_FetalKidney/Human_FetalKidney.RData")
