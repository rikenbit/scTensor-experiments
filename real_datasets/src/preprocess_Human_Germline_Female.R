source("src/Functions.R")

########################################################
# 1. Gene expression matrix
########################################################

# Data Loading
GSE86146 <- read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306011_F_5W_embryo1_and_2_gene_expression.txt",
	header=TRUE)
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306011_F_5W_embryo1_and_2_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306012_F_7W_embryo1_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306013_F_8W_embryo1_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2295850_F_10W_embryo1_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306014_F_12W_embryo1_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2295851_F_14W_embryo1_1_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2295852_F_14W_embryo1_2_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2295853_F_14W_embryo1_3_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2295854_F_18W_embryo1_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306015_F_18W_embryo2_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306016_F_20W_embryo1_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306017_F_20W_embryo2_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306018_F_23W_embryo1_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306019_F_23W_embryo2_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306020_F_24W_embryo1_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306021_F_24W_embryo2_gene_expression.txt",
	header=TRUE))
GSE86146 <- cbind(GSE86146, read.table("data/Human_Germline_Female/GSE86146_RAW/GSM2306022_F_26W_embryo1_gene_expression.txt",
	header=TRUE))

rownames.GSE86146 <- as.character(GSE86146[,1])
GSE86146 <- GSE86146[, setdiff(1:ncol(GSE86146),
	grep("Gene", colnames(GSE86146)))]
rownames(GSE86146) <- rownames.GSE86146

# Gene ID conversion
GSE86146_Female <- GSE86146
LefttoRight <- select(Homo.sapiens,
  column=c("SYMBOL", "ENTREZID"),
  keytype="SYMBOL",
  keys=rownames(GSE86146_Female))

GSE86146_Female <- convertToNCBIGeneID(
	GSE86146_Female,
	rownames(GSE86146_Female),
	LefttoRight)

# cell type label
label.GSE86146 <- read.table("data/Human_Germline_Female/Cluster_Germline.txt", header=TRUE)
label.GSE86146_Female <- intersect(colnames(GSE86146_Female), label.GSE86146[, 1])

GSE86146_Female <- GSE86146_Female[, label.GSE86146_Female]

names(label.GSE86146_Female) <- as.character(label.GSE86146[sapply(label.GSE86146_Female, function(x){which(label.GSE86146[, 1] == x)}), 2])

label.GSE86146_Female[which(names(label.GSE86146_Female) == "Female_FGC_1")] <- brewer.pal(11, "Spectral")[1]
label.GSE86146_Female[which(names(label.GSE86146_Female) == "Female_FGC_2")] <- brewer.pal(11, "Spectral")[2]
label.GSE86146_Female[which(names(label.GSE86146_Female) == "Female_FGC_3")] <- brewer.pal(11, "Spectral")[3]
label.GSE86146_Female[which(names(label.GSE86146_Female) == "Female_FGC_4")] <- brewer.pal(11, "Spectral")[4]
label.GSE86146_Female[which(names(label.GSE86146_Female) == "Female_Soma_1")] <- brewer.pal(11, "Spectral")[8]
label.GSE86146_Female[which(names(label.GSE86146_Female) == "Female_Soma_2")] <- brewer.pal(11, "Spectral")[9]
label.GSE86146_Female[which(names(label.GSE86146_Female) == "Female_Soma_3")] <- brewer.pal(11, "Spectral")[10]
label.GSE86146_Female[which(names(label.GSE86146_Female) == "Female_Soma_4")] <- brewer.pal(11, "Spectral")[11]

hvg.stats <- HVG(GSE86146_Female)
target <- which(hvg.stats$pval < 1E-7)
GSE86146_Female.HVG <- GSE86146_Female[target, ]

input <- GSE86146_Female.HVG

########################################################
# 2. Celltype label vector
########################################################
celltypes <- label.GSE86146_Female

########################################################
# 3. Ligand-Receptor pairs used in the original paper
########################################################
db <- c(
	"CELLPHONEDB", "SINGLECELLSIGNALR", "DLRP", "IUPHAR", "HPMR")

tmp1 <- select(LRBase.Hsa.eg.db, columns=c("GENEID_L", "GENEID_R", "SOURCEDB"), keytype="GENEID_L", keys=rownames(input))
tmp2 <- select(LRBase.Hsa.eg.db, columns=c("GENEID_L", "GENEID_R", "SOURCEDB"), keytype="GENEID_R", keys=rownames(input))

tmp1 <- tmp1[which(tmp1$SOURCEDB %in% db), c("GENEID_L", "GENEID_R")]
tmp2 <- tmp2[which(tmp2$SOURCEDB %in% db), c("GENEID_L", "GENEID_R")]

LR <- unique(rbind(tmp1, tmp2))

########################################################
# Saving
########################################################
save(input, LR, celltypes,
	file="data/Human_Germline_Female/Human_Germline_Female.RData")
