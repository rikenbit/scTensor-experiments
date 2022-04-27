source("src/Functions.R")

# Download
td <- tempdir()
destfile <- paste0(td, "E-MTAB-6583.processed.1.zip")
download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6583/E-MTAB-6583.processed.1.zip",
	destfile=destfile)
unzip(destfile, exdir=td)

# 問い合わせ中
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6686/samples/

# CTRL(Lgr5), STATE1(Lgr5), CTRL(Lgr6), STATE1(Lgr6),
# CTRL (Bulk),STATE1(Bulk)
# からSTATE1(Lgr5)/STATE1(Lgr6)/24h(Bulk)のCCIを見る解析

# Data Loading
E_MTAB_6583 <- read.delim(paste0(td, "/201603171136_v2.6_seq_WND.txt"),
	sep="\t", header=TRUE, stringsAsFactor=FALSE)
rownames(E_MTAB_6583) <- E_MTAB_6583[,1]

# Label
label.E_MTAB_6583_Lgr5 <- read.delim("data/Mouse_WoundHealing/201603171136_v2.6_cl_Lgr5_wnd_ctrl_sel.txt",
	header=FALSE, stringsAsFactor=FALSE)
label.E_MTAB_6583_Lgr6 <- read.delim("data/Mouse_WoundHealing/201603171136_v2.6_cl_Lgr6_wnd_ctrl_sel.txt",
	header=FALSE, stringsAsFactor=FALSE)

# Preprocess
E_MTAB_6583_Lgr5 <- E_MTAB_6583[,
	gsub("-", ".", paste0("X", label.E_MTAB_6583_Lgr5[,1]))]
E_MTAB_6583_Lgr6 <- E_MTAB_6583[,
	gsub("-", ".", paste0("X", label.E_MTAB_6583_Lgr6[,1]))]

# Label Color Setting
label.E_MTAB_6583_Lgr5 <- label.E_MTAB_6583_Lgr5[,2]
label.E_MTAB_6583_Lgr5[which(label.E_MTAB_6583_Lgr5 == 3)] <- "State 0"
label.E_MTAB_6583_Lgr5[which(label.E_MTAB_6583_Lgr5 == 2)] <- "State 1A"
label.E_MTAB_6583_Lgr5[which(label.E_MTAB_6583_Lgr5 == 0)] <- "State 1B"
label.E_MTAB_6583_Lgr5[which(label.E_MTAB_6583_Lgr5 == 1)] <- "State 2A"
label.E_MTAB_6583_Lgr5[which(label.E_MTAB_6583_Lgr5 == 4)] <- "State 2B"

label.E_MTAB_6583_Lgr6 <- label.E_MTAB_6583_Lgr6[,2]
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == 5)] <- "State 1B"
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == 1)] <- "State 1A"
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == 4)] <- "State 2"
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == 2)] <- "State 3A"
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == 0)] <- "State 3B"
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == 3)] <- "uHF-like"

names(label.E_MTAB_6583_Lgr5) <- label.E_MTAB_6583_Lgr5
names(label.E_MTAB_6583_Lgr6) <- label.E_MTAB_6583_Lgr6

# ctrl : 黄色
label.E_MTAB_6583_Lgr5[which(label.E_MTAB_6583_Lgr5 == "ctrl")] <-
	brewer.pal(8, "Set2")[6]
# State0 : オレンジ
label.E_MTAB_6583_Lgr5[which(label.E_MTAB_6583_Lgr5 == "State 0")] <-
	brewer.pal(12, "Paired")[8]
# State1A : 薄紫
label.E_MTAB_6583_Lgr5[which(label.E_MTAB_6583_Lgr5 == "State 1A")] <-
	brewer.pal(12, "Paired")[9]
# State1B : 紫
label.E_MTAB_6583_Lgr5[which(label.E_MTAB_6583_Lgr5 == "State 1B")] <-
	brewer.pal(12, "Paired")[10]
# State2A : 薄緑
label.E_MTAB_6583_Lgr5[which(label.E_MTAB_6583_Lgr5 == "State 2A")] <-
	brewer.pal(12, "Paired")[3]
# State2B : 緑
label.E_MTAB_6583_Lgr5[which(label.E_MTAB_6583_Lgr5 == "State 2B")] <-
	brewer.pal(12, "Paired")[4]

# ctrl - IFE : 黄色1
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == "ctrl - IFE")] <-
	brewer.pal(8, "Set2")[6]
# ctrl - IST : 黄色2
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == "ctrl - IST")] <-
	brewer.pal(9, "Set1")[6]
# State 1A : 薄紫
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == "State 1A")] <-
	brewer.pal(12, "Paired")[9]
# State 1B : 紫
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == "State 1B")] <-
	brewer.pal(12, "Paired")[10]
# State 2 : 緑
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == "State 2")] <-
	brewer.pal(12, "Paired")[4]
# State 3A : 薄赤
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == "State 3A")] <-
	brewer.pal(12, "Paired")[5]
# State 3B : 赤
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == "State 3B")] <-
	brewer.pal(12, "Paired")[6]
# uHF-like : 茶色
label.E_MTAB_6583_Lgr6[which(label.E_MTAB_6583_Lgr6 == "uHF-like")] <-
	brewer.pal(12, "Paired")[12]

# t-SNE corrdinates
res.tsne.E_MTAB_6583_Lgr5 <- read.delim("data/Mouse_WoundHealing/201603171136_v2.6_tsne_Lgr5_wnd_ctrl_sel.txt",
	sep="\t", header=TRUE, stringsAsFactor=FALSE)
res.tsne.E_MTAB_6583_Lgr6 <- read.delim("data/Mouse_WoundHealing/201603171136_v2.6_tsne_Lgr6_wnd_ctrl_sel.txt",
	sep="\t", header=TRUE, stringsAsFactor=FALSE)
rownames(res.tsne.E_MTAB_6583_Lgr5) <- res.tsne.E_MTAB_6583_Lgr5[,1]
rownames(res.tsne.E_MTAB_6583_Lgr6) <- res.tsne.E_MTAB_6583_Lgr6[,1]
res.tsne.E_MTAB_6583_Lgr5 <- res.tsne.E_MTAB_6583_Lgr5[,c("x", "y")]
res.tsne.E_MTAB_6583_Lgr6 <- res.tsne.E_MTAB_6583_Lgr6[,c("x", "y")]





# HVG
hvg.stats <- HVG(cbind(E_MTAB_6583_Lgr5, E_MTAB_6583_Lgr6))
target <- which(hvg.stats$pval < 1E-1)

# hist(hvg.stats$pval)

# plot(log10(hvg.stats$m+1), log10(hvg.stats$cv2+1), xlim=quantile(log10(hvg.stats$m+1))[c(1,5)], ylim=quantile(log10(hvg.stats$cv2+1))[c(1,5)])
# par(new=TRUE)
# plot(log10(hvg.stats$m[target]+1), log10(hvg.stats$cv2[target]+1), col="red", xlim=quantile(log10(hvg.stats$m+1))[c(1,5)], ylim=quantile(log10(hvg.stats$cv2+1))[c(1,5)])

E_MTAB_6583_Lgr5.HVG <- E_MTAB_6583_Lgr5[target, ]
E_MTAB_6583_Lgr6.HVG <- E_MTAB_6583_Lgr6[target, ]





# PCA
res.pca.E_MTAB_6583_Lgr5 <- prcomp(log10(E_MTAB_6583_Lgr5+1))
res.pca.E_MTAB_6583_Lgr6 <- prcomp(log10(E_MTAB_6583_Lgr6+1))
res.pca.E_MTAB_6583_Lgr5.HVG <- prcomp(log10(E_MTAB_6583_Lgr5.HVG+1))
res.pca.E_MTAB_6583_Lgr6.HVG <- prcomp(log10(E_MTAB_6583_Lgr6.HVG+1))


# ID Conversion: Gene Symbol => NCBI Gene ID
LtoR <- Mus.musculus
E_MTAB_6583_Lgr5 <- scTGIF::convertRowID(, E_MTAB_6583_Lgr5)
E_MTAB_6583_Lgr6 <- scTGIF::convertRowID(Mus.musculus, E_MTAB_6583_Lgr6)




# 保存
save(E_MTAB_6583_Lgr5, E_MTAB_6583_Lgr5.HVG,
	label.E_MTAB_6583_Lgr5,
  res.pca.E_MTAB_6583_Lgr5, res.pca.E_MTAB_6583_Lgr5.HVG,
  res.tsne.E_MTAB_6583_Lgr5,
  res.umap.E_MTAB_6583_Lgr5, res.umap.E_MTAB_6583_Lgr5.HVG,
  file="../Data/Mouse/WoundHealing_Lgr5/WoundHealing_Lgr5.RData")

save(E_MTAB_6583_Lgr6, E_MTAB_6583_Lgr6.HVG,
	label.E_MTAB_6583_Lgr6,
  res.pca.E_MTAB_6583_Lgr6, res.pca.E_MTAB_6583_Lgr6.HVG,
  res.tsne.E_MTAB_6583_Lgr6,
  res.umap.E_MTAB_6583_Lgr6, res.umap.E_MTAB_6583_Lgr6.HVG,
  file="../Data/Mouse/WoundHealing_Lgr6/WoundHealing_Lgr6.RData")


















# # プロット（Lgr5）
# png(file="plot/WoundHealing_Lgr5/tSNE.png", width=750, height=750)
# plot(res.tsne.E_MTAB_6583_Lgr5, col=label.E_MTAB_6583_Lgr5, pch=16, xlab="Dim1", ylab="Dim2", cex=2)
# dev.off()

# png(file="plot/WoundHealing_Lgr5/UMAP.png", width=750, height=750)
# plot(res.umap.E_MTAB_6583_Lgr5$layout, col=label.E_MTAB_6583_Lgr5, pch=16, xlab="Dim1", ylab="Dim2", cex=2)
# dev.off()

# png(file="plot/WoundHealing_Lgr5/UMAP_HVG.png", width=750, height=750)
# plot(res.umap.E_MTAB_6583_Lgr5.HVG$layout, col=label.E_MTAB_6583_Lgr5, pch=16, xlab="Dim1", ylab="Dim2", cex=2)
# dev.off()

# # プロット（Lgr6）
# png(file="plot/WoundHealing_Lgr6/tSNE.png", width=750, height=750)
# plot(res.tsne.E_MTAB_6583_Lgr6, col=label.E_MTAB_6583_Lgr6, pch=16, xlab="Dim1", ylab="Dim2", cex=2)
# dev.off()

# png(file="plot/WoundHealing_Lgr6/UMAP.png", width=750, height=750)
# plot(res.umap.E_MTAB_6583_Lgr6$layout, col=label.E_MTAB_6583_Lgr6, pch=16, xlab="Dim1", ylab="Dim2", cex=2)
# dev.off()

# png(file="plot/WoundHealing_Lgr6/UMAP_HVG.png", width=750, height=750)
# plot(res.umap.E_MTAB_6583_Lgr6.HVG$layout, col=label.E_MTAB_6583_Lgr6, pch=16, xlab="Dim1", ylab="Dim2", cex=2)
# dev.off()
