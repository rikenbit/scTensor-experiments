source("src/Functions.R")

# Download
td <- tempdir()
destfile1 <- paste0(td, "GSE102827_merged_all_raw.csv.gz")
destfile2 <- paste0(td, "GSE102827_cell_type_assignments.csv.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102827/suppl/GSE102827%5Fmerged%5Fall%5Fraw%2Ecsv%2Egz",
	destfile=destfile1)
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102827/suppl/GSE102827%5Fcell%5Ftype%5Fassignments%2Ecsv%2Egz",
	destfile=destfile2)

# Data Loading
GSE102827 <- read.delim(destfile1, sep=",", header=TRUE)
label.GSE102827 <- read.delim(destfile2, sep=",", header=TRUE,
	stringsAsFactor=FALSE)

row.names <- GSE102827[,1]
GSE102827 <- GSE102827[, 2:ncol(GSE102827)]
rownames(GSE102827) <- row.names

# Re-arrange
order.fig2 <- c(
	"ExcL23", "ExcL4", "ExcL5_1", "ExcL5_2", "ExcL5_3", "ExcL6",
	"Int_Pv", "Int_Sst_1", "Int_Sst_2", "Int_Vip", "Int_Npy", "Int_Cck",
	"Olig_1", "Olig_2", "Olig_3", "Olig_4", "Olig_5", "Olig_6", "Olig_7",
	"OPC_1", "OPC_2",
	"Endo_1", "Endo_2",
	"SM_1", "SM_2",
	"Micro_1", "Micro_2",
	"Astro", "Pericyte", "Macrophage")
target <- intersect(
	unique(colnames(GSE102827)),
	unique(label.GSE102827$X[which(label.GSE102827$celltype %in% order.fig2)]))
GSE102827 <- GSE102827[, target]
target2 <- unlist(sapply(target, function(x){
	which(label.GSE102827$X == x)}))
label.GSE102827 <- label.GSE102827$celltype[target2]

# Label Color
names(label.GSE102827) <- label.GSE102827

# Excitatory : 青
label.GSE102827[
	which(label.GSE102827 == "ExcL23")] <- "#5E4FA2"
label.GSE102827[
	which(label.GSE102827 == "ExcL4")] <- "#5E4FA3"
label.GSE102827[
	which(label.GSE102827 == "ExcL5_1")] <- "#5E4FA4"
label.GSE102827[
	which(label.GSE102827 == "ExcL5_2")] <- "#5E4FA5"
label.GSE102827[
	which(label.GSE102827 == "ExcL5_3")] <- "#5E4FA6"
label.GSE102827[
	which(label.GSE102827 == "ExcL6")] <- "#5E4FA7"

# Interneurons : ピンク
label.GSE102827[
	which(label.GSE102827 == "Int_Pv")] <- "#D53E4F"
label.GSE102827[
	which(label.GSE102827 == "Int_Sst_1")] <- "#D53E4E"
label.GSE102827[
	which(label.GSE102827 == "Int_Sst_2")] <- "#D53E4D"
label.GSE102827[
	which(label.GSE102827 == "Int_Vip")] <- "#D53E4C"
label.GSE102827[
	which(label.GSE102827 == "Int_Npy")] <- "#D53E4B"
label.GSE102827[
	which(label.GSE102827 == "Int_Cck")] <- "#D53E4A"

# Oligodendrocytes : 青色
label.GSE102827[
	which(label.GSE102827 == "Olig_1")] <- "#3288BD"
label.GSE102827[
	which(label.GSE102827 == "Olig_2")] <- "#3288BE"
label.GSE102827[
	which(label.GSE102827 == "Olig_3")] <- "#3288BF"
label.GSE102827[
	which(label.GSE102827 == "Olig_4")] <- "#3288BC"
label.GSE102827[
	which(label.GSE102827 == "Olig_5")] <- "#3288BB"
label.GSE102827[
	which(label.GSE102827 == "Olig_6")] <- "#3288BA"
label.GSE102827[
	which(label.GSE102827 == "Olig_7")] <- "#3288CA"
label.GSE102827[
	which(label.GSE102827 == "OPC_1")] <- "#3288CB"
label.GSE102827[
	which(label.GSE102827 == "OPC_2")] <- "#3288CC"

# Endothelial_SmoothMuscle : 青緑
label.GSE102827[
	which(label.GSE102827 == "Endo_1")] <- "#66C2A5"
label.GSE102827[
	which(label.GSE102827 == "Endo_2")] <- "#66C2A4"
label.GSE102827[
	which(label.GSE102827 == "SM_1")] <- "#66C2A3"
label.GSE102827[
	which(label.GSE102827 == "SM_2")] <- "#66C2A2"

# Microglia : オレンジ
label.GSE102827[
	which(label.GSE102827 == "Micro_1")] <- "#F46D43"
label.GSE102827[
	which(label.GSE102827 == "Micro_2")] <- "#F46D42"

# Astrocytes : 山吹色
label.GSE102827[
	which(label.GSE102827 == "Astro")] <- "#FDAE61"

# Mural = Pericytes? : フカビドリ
label.GSE102827[
	which(label.GSE102827 == "Pericyte")] <- "#006837"

# Macrophages : 茶色
# Macrophage"
label.GSE102827[
	which(label.GSE102827 == "Macrophage")] <- "#543005"

# Gene IDに書き換え（計算待ち）
LtoR <- select(Mus.musculus,
	columns=c("SYMBOL", "GENEID"),
	keys=rownames(GSE102827),
	keytype="SYMBOL")
obj.scTGIF <- convertRowID(GSE102827, rownames(GSE102827), LtoR)
GSE102827 <- obj.scTGIF$output

# HVGに絞る
hvg.stats <- HVG(GSE102827)
GSE102827.HVG <- GSE102827[which(hvg.stats$pval < 1E-1), ]

# hist(hvg.stats$pval)

# plot(log10(hvg.stats$m+1), log10(hvg.stats$cv2+1), xlim=quantile(log10(hvg.stats$m+1))[c(1,5)], ylim=quantile(log10(hvg.stats$cv2+1))[c(1,5)])
# par(new=TRUE)
# plot(log10(hvg.stats$m[target]+1), log10(hvg.stats$cv2[target]+1), col="red", xlim=quantile(log10(hvg.stats$m+1))[c(1,5)], ylim=quantile(log10(hvg.stats$cv2+1))[c(1,5)])

res.pca.GSE102827.HVG <- prcomp_irlba(log10(GSE102827.HVG + 1), n=30)
res.umap.GSE102827.HVG <- umap(t(log10(GSE102827.HVG + 1)))

# Ligand-Receptor pairs
db <- c(
	"ENSEMBL_DLRP", "ENSEMBL_IUPHAR", "ENSEMBL_HPMR",
	"NCBI_DLRP", "NCBI_IUPHAR", "NCBI_HPMR")
setAnnotationHubOption("CACHE", getwd())
ah <- AnnotationHub()
LRBase.Mmu.eg.db <- query(ah,
	c("LRBaseDb", "Mus musculus", "v002"))[[1]]
LRBase.Mmu.eg.db <- LRBaseDb(LRBase.Mmu.eg.db)

tmp1 <- select(LRBase.Mmu.eg.db,
	columns=c("GENEID_L", "GENEID_R", "SOURCEDB"),
	keytype="GENEID_L", keys=rownames(GSE102827.HVG))
tmp2 <- select(LRBase.Mmu.eg.db,
	columns=c("GENEID_L", "GENEID_R", "SOURCEDB"),
	keytype="GENEID_R", keys=rownames(GSE102827.HVG))

tmp1 <- tmp1[which(tmp1$SOURCEDB %in% db),
	c("GENEID_L", "GENEID_R")]
tmp2 <- tmp2[which(tmp2$SOURCEDB %in% db),
	c("GENEID_L", "GENEID_R")]

LR <- unique(rbind(tmp1, tmp2))

# 保存
input = GSE102827.HVG
celltypes = label.GSE102827
save(input, LR, celltypes,
	file="data/Mouse_VisualCortex/Mouse_VisualCortex.RData")

# # プロット
# plot(res.umap.GSE102827.HVG$layout, pch=16, xlab="Dim1", ylab="Dim2", cex=0.5)

# plot(res.umap.GSE102827.HVG$layout, pch=16, xlab="Dim1", ylab="Dim2", cex=0.5, col=label.GSE102827)
