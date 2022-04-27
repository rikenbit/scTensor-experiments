source("src/Functions.R")

# Download
destfile <- paste0(tempdir(), "GSE103322_HNSCC_all_data.txt.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/GSE103322%5FHNSCC%5Fall%5Fdata%2Etxt%2Egz", destfile=destfile)

# Preprocess
GSE103322 <- read.delim(destfile, sep="\t",
	stringsAsFactor=FALSE, row.names=1)

malignants <- c("HN25", "HN26", "HN28", "HNSCC",
  "HNSCC12", "HNSCC13", "HNSCC16", "HNSCC17", "HNSCC18",
  "HNSCC20", "HNSCC22", "HNSCC24", "HNSCC25", "HNSCC26",
  "HNSCC28", "HNSCC5", "HNSCC6", "HNSCC7")
malignant_patterns <- paste0("^", malignants, "_")
malignants_colors <- c(
  brewer.pal(11, "Set3"),
  brewer.pal(9, "Pastel1"))

label.GSE103322 <- as.character(GSE103322[3,])
label.celltype <- as.character(GSE103322[5,])
for(i in seq(malignants)){
  label.GSE103322[
  intersect(
    which(label.GSE103322 == 1),
    grep(malignant_patterns[i], colnames(GSE103322)))] <- malignants[i]

}
label.GSE103322[which(label.celltype == "-Fibroblast")] <- "CAF"
label.GSE103322[which(label.celltype == "B cell")] <- "B cell"
label.GSE103322[which(label.celltype == "Dendritic")] <- "Dendritic"
label.GSE103322[which(label.celltype == "Endothelial")] <- "Endothelial"
label.GSE103322[which(label.celltype == "Fibroblast")] <- "Fibroblast"
label.GSE103322[which(label.celltype == "Macrophage")] <- "Macrophage"
label.GSE103322[which(label.celltype == "Mast")] <- "Mast"
label.GSE103322[which(label.celltype == "myocyte")] <- "myocyte"
label.GSE103322[which(label.celltype == "T cell")] <- "T cell"
names(label.GSE103322) <- label.GSE103322

target <- which(label.GSE103322 %ni% c(0,1))
label.GSE103322 <- label.GSE103322[target]
GSE103322 <- as.matrix(GSE103322[6:ncol(GSE103322), target])
rn <- gsub("'", "", rownames(GSE103322))
cn <- colnames(GSE103322)
newGSE103322 <- matrix(0, nrow=length(rn), ncol=length(cn))
for(i in 1:ncol(GSE103322)){
  newGSE103322[,i] <- as.numeric(GSE103322[,i])
}
rownames(newGSE103322) <- rn
colnames(newGSE103322) <- cn
GSE103322 <- newGSE103322

# 色情報
label.GSE103322[which(label.GSE103322 == "T cell")] <- rgb(35/255,68/255,147/255)
label.GSE103322[which(label.GSE103322 == "B cell")] <- rgb(30/255,107/255,180/255)
label.GSE103322[which(label.GSE103322 == "Macrophage")] <- rgb(0/255,133/255,195/255)
label.GSE103322[which(label.GSE103322 == "Mast")] <- rgb(0/255,214/255,239/255)
label.GSE103322[which(label.GSE103322 == "Dendritic")] <- rgb(99/255,211/255,172/255)
label.GSE103322[which(label.GSE103322 == "Endothelial")] <- rgb(202/255,220/255,63/255)
label.GSE103322[which(label.GSE103322 == "CAF")] <- rgb(175/255,15/255,21/255)
label.GSE103322[which(label.GSE103322 == "Fibroblast")] <- rgb(255/255,7/255,32/255)
label.GSE103322[which(label.GSE103322 == "myocyte")] <- rgb(255/255,111/255,28/255)
for(i in seq(malignants)){
  label.GSE103322[
    which(label.GSE103322 == malignants[i])] <- malignants_colors[i]
}

# HVG
hvg.stats <- HVG(GSE103322)
target <- which(hvg.stats$pval < 1E-1)
GSE103322.HVG <- GSE103322[target, ]

# PCA
res.pca.GSE103322 <- prcomp(log10(GSE103322 + 1), scale=TRUE)
res.pca.GSE103322.HVG <- prcomp(log10(GSE103322.HVG + 1), scale=TRUE)

# UMAP
set.seed(123)
res.umap.GSE103322 <- umap(t(log10(GSE103322+1)))
set.seed(123)
res.umap.GSE103322.HVG <- umap(t(log10(GSE103322.HVG+1)))

# Gene ID conversion
LefttoRight <- select(Homo.sapiens,
  column=c("SYMBOL", "ENTREZID"),
  keytype="SYMBOL",
  keys=rownames(GSE103322))

GSE103322 <- convertRowID(
  GSE103322,
  rownames(GSE103322),
  LefttoRight)$output

db <- c("DLRP", "IUPHAR", "HPMR")
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
input <- GSE103322
celltypes <- label.GSE103322
save(input, LR, celltypes,
  file="data/Human_HeadandNeckCancer/Human_HeadandNeckCancer.RData")
  