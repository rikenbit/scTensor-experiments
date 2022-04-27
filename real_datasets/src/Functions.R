library("igraph")
library("scTensor")
library("scTGIF")
library("scales")
library("RColorBrewer")
library("ROCR")
library("nnTensor")
library("rTensor")
library("Matrix")
library("Homo.sapiens")
library("Mus.musculus")
library("DESeq2")
library("statmod")
# library("LRBase.Hsa.eg.db")
# library("LRBase.Mmu.eg.db")
library("LRBaseDbi")
library("AnnotationHub")
library("ggplot2")
library("viridis")
library("dbscan")
library("mclust")
library("uwot")
library("Rtsne")
library("irlba")

Methods <- c("labelpermutation_tensor",
"labelpermutation",
"labelpermutation2_tensor",
"labelpermutation2",
"halpern_tensor",
"halpern",
"cabelloaguilar_tensor",
"cabelloaguilar",
"previous_sctensor",
"sctensor")

BinMethods <- c(
"labelpermutation",
"labelpermutation2",
"halpern",
"cabelloaguilar",
"previous_sctensor",
"sctensor")

Samples <- c("Human_FetalKidney",
    "Human_NicotinehESCs_Nicotine",
    "Human_Germline_Female",
    "Human_HeadandNeckCancer",
    "Mouse_Uterus",
    "Mouse_VisualCortex")

Label.Samples <- c("Human\nFetalKidney",
    "Human\nNicotinehESCs\nNicotine",
    "Human\nGermline\nFemale",
    "Human\nHeadandNeckCancer",
    "Mouse\nUterus",
    "Mouse\nVisualCortex")


'%ni%' <- Negate('%in%')

options(timeout=1e10)
LogCPM <- function(input){
    libsize <- colSums(input)
    cpm <- median(libsize) * t(t(input) / libsize)
    log10(cpm + 1)
}

# Check 2.4.1
if(packageVersion("scTensor") != "2.4.1"){
stop("scTensor is not 2.4.1!")
}

# Plot Ground Truth
mylayout <- function(l){
    rbind(
        cbind(2, rev(seq(l))),
        cbind(1, rev(seq(l))))
}

bipartiteGraph <- function(cif, x){
    # Setting
    eachcif <- cif[[x]]
    l <- length(eachcif$LPattern)
    inc <- outer(eachcif$RPattern, eachcif$LPattern )
    rownames(inc) <- paste0("Receptor-", seq(l))
    colnames(inc) <- paste0("Ligand-", seq(l))
    g <- graph_from_incidence_matrix(inc)
    mynodecolor <- c(rgb(0,0,1,0.5), rgb(1,0,0,0.5))
    myedgecolor <- c(
        brewer.pal(8, "Dark2"),
        brewer.pal(9, "Set1"),
        brewer.pal(7, "Pastel1")
        )
    plot(g, layout = mylayout(l),
        vertex.color=mynodecolor[V(g)$type + 1],
        edge.color=myedgecolor[x-1])
}

mypanel <- function(l){
    if((1 < l) && (l <= 5)){
        layout(t(seq(l)))
    }
    if((6 < l) && (l <= 10)){
        layout(rbind(t(1:5), t(6:10)))
    }
    if((11 < l) && (l <= 15)){
        layout(rbind(t(1:5), t(6:10), t(11:15)))
    }
    if((16 < l) && (l <= 20)){
        layout(rbind(t(1:5), t(6:10), t(11:15), t(16:20)))
    }
    if((21 < l) && (l <= 25)){
        layout(rbind(t(1:5), t(6:10), t(11:15), t(16:20), t(21:25)))
    }
}

mywidth <- function(l){
    rep(seq(5), 5)[l] * 300
}

myheight <- function(l){
    if((1 < l) && (l <= 5)){
        height <- 300
    }
    if((6 < l) && (l <= 10)){
        height <- 2 * 300
    }
    if((11 < l) && (l <= 15)){
        height <- 3 * 300
    }
    if((16 < l) && (l <= 20)){
        height <- 4 * 300
    }
    if((21 < l) && (l <= 25)){
        height <- 5 * 300
    }
    height
}

# For ROC, AUC, Binarization, F-measure
tensor.methods <- c(
    "labelpermutation_tensor", "labelpermutation2_tensor",
    "halpern_tensor", "cabelloaguilar_tensor")
pval.methods <- c(
    "labelpermutation", "labelpermutation2",
    "halpern", "cabelloaguilar")

# Memory
Plot_Memory <- function(){
    df <- aggregateMemory()
    # Plot
    gg <- ggplot(df, aes(x=Methods, y=Samples, fill=GB))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(Samples~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="GB", limits=c(0, max(df$GB)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=1, size=9))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave("plot/Memory.png", gg, dpi=120, width=10, height=12)
}

aggregateMemory <- function(){
    out <- c()
    nM <- length(Methods)
    nS <- length(Samples)
    for(i in 1:nM){
        for(j in 1:nS){
            if(Methods[i] %in% c("sctensor", "previous_sctensor")){
                inputfile <- paste0("benchmarks/",
                    Methods[i], "_1_",
                    Samples[j], ".txt")
            }else{
                inputfile <- paste0("benchmarks/",
                    gsub("_tensor", "", Methods[i]), "_",
                    Samples[j], ".txt")
            }
            memory <- read.delim(inputfile)
            out <- c(out, memory$max_rss)
        }
    }
    tmp <- as.vector(
            sapply(Methods, function(x){
                rep(x, length=nS)}))
    df <- data.frame(
        Methods=tmp,
        Samples=Label.Samples,
        GB=out)
    df$Methods <- factor(df$Methods,
        levels=unique(df$Methods))
    df$Samples <- factor(df$Samples,
        levels=unique(df$Samples))
    df
}

# Time
Plot_Time <- function(){
    df <- aggregateTime()
    # Plot
    gg <- ggplot(df, aes(x=Methods, y=Samples, fill=Time))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(Samples~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="Hour", limits=c(0, max(df$Time)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=1, size=9))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave("plot/Time.png", gg, dpi=120, width=10, height=12)
}

aggregateTime <- function(){
    out <- c()
    nM <- length(Methods)
    nS <- length(Samples)
    for(i in 1:nM){
        for(j in 1:nS){
            if(Methods[i] %in% c("sctensor", "previous_sctensor")){
                inputfile <- paste0("benchmarks/",
                    Methods[i], "_1_",
                    Samples[j], ".txt")
            }else{
                inputfile <- paste0("benchmarks/",
                    gsub("_tensor", "", Methods[i]), "_",
                    Samples[j], ".txt")
            }
            time <- read.delim(inputfile)
            out <- c(out, time$s / 60 / 60)
        }
    }
    tmp <- as.vector(
            sapply(Methods, function(x){
                rep(x, length=nS)}))
    df <- data.frame(
        Methods=tmp,
        Samples=Label.Samples,
        Time=out)
    df$Methods <- factor(df$Methods,
        levels=unique(df$Methods))
    df$Samples <- factor(df$Samples,
        levels=unique(df$Samples))
    df
}

# F
Plot_F <- function(){
    df <- aggregateF()
    # Plot
    gg <- ggplot(df, aes(x=BinMethods, y=Samples, fill=F))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(Samples~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="F-measure", limits=c(0, max(df$F)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=1, size=9))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave("plot/F.png", gg, dpi=120, width=10, height=12)
}

aggregateF <- function(){
    out <- c()
    nM <- length(BinMethods)
    nS <- length(Samples)
    for(i in 1:nM){
        for(j in 1:nS){
            inputfile <- paste0("output/",
                BinMethods[i], "/F/", Samples[j], ".RData")
            load(inputfile)
            f <- unlist(f)
            target <- which(!is.nan(f))
            out = c(out, mean(f[target]))
        }
    }
    tmp <- as.vector(
            sapply(BinMethods, function(x){
                rep(x, length=nS)}))
    df <- data.frame(
        BinMethods=tmp,
        Samples=Label.Samples,
        F=out)
    df$BinMethods <- factor(df$BinMethods,
        levels=unique(df$BinMethods))
    df$Samples <- factor(df$Samples,
        levels=unique(df$Samples))
    df
}

# MCC
Plot_MCC <- function(){
    df <- aggregateMCC()
    # Plot
    gg <- ggplot(df, aes(x=BinMethods, y=Samples, fill=MCC))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(Samples~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="MCC", limits=c(0, max(df$MCC)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=1, size=9))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave("plot/MCC.png", gg, dpi=120, width=10, height=12)
}

aggregateMCC <- function(){
    out <- c()
    nM <- length(BinMethods)
    nS <- length(Samples)
    for(i in 1:nM){
        for(j in 1:nS){
            inputfile <- paste0("output/",
                BinMethods[i], "/MCC/", Samples[j], ".RData")
            load(inputfile)
            mcc <- unlist(mcc)
            target <- which(!is.nan(mcc))
            out = c(out, mean(mcc[target]))
        }
    }
    tmp <- as.vector(
            sapply(BinMethods, function(x){
                rep(x, length=nS)}))
    df <- data.frame(
        BinMethods=tmp,
        Samples=Label.Samples,
        MCC=out)
    df$BinMethods <- factor(df$BinMethods,
        levels=unique(df$BinMethods))
    df$Samples <- factor(df$Samples,
        levels=unique(df$Samples))
    df$MCC[which(df$MCC < 0)] <- 0
    df
}

# FPR
Plot_FPR <- function(){
    df <- aggregateFPR()
    # Plot
    gg <- ggplot(df, aes(x=BinMethods, y=Samples, fill=FPR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(Samples~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="FPR", limits=c(0, max(df$FPR)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=1, size=9))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave("plot/FPR.png", gg, dpi=120, width=10, height=12)
}

aggregateFPR <- function(){
    out <- c()
    nM <- length(BinMethods)
    nS <- length(Samples)
    for(i in 1:nM){
        for(j in 1:nS){
            inputfile <- paste0("output/",
                BinMethods[i], "/FPR/", Samples[j], ".RData")
            load(inputfile)
            fpr <- unlist(fpr)
            target <- which(!is.nan(fpr))
            out = c(out, mean(fpr[target]))
        }
    }
    tmp <- as.vector(
            sapply(BinMethods, function(x){
                rep(x, length=nS)}))
    df <- data.frame(
        BinMethods=tmp,
        Samples=Label.Samples,
        FPR=out)
    df$BinMethods <- factor(df$BinMethods,
        levels=unique(df$BinMethods))
    df$Samples <- factor(df$Samples,
        levels=unique(df$Samples))
    df$FPR[which(df$FPR < 0)] <- 0
    df
}

# FNR
Plot_FNR <- function(){
    df <- aggregateFNR()
    # Plot
    gg <- ggplot(df, aes(x=BinMethods, y=Samples, fill=FNR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(Samples~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="FNR", limits=c(0, max(df$FNR)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=1, size=9))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave("plot/FNR.png", gg, dpi=120, width=10, height=12)
}

aggregateFNR <- function(){
    out <- c()
    nM <- length(BinMethods)
    nS <- length(Samples)
    for(i in 1:nM){
        for(j in 1:nS){
            inputfile <- paste0("output/",
                BinMethods[i], "/FNR/", Samples[j], ".RData")
            load(inputfile)
            fnr <- unlist(fnr)
            target <- which(!is.nan(fnr))
            out = c(out, mean(fnr[target]))
        }
    }
    tmp <- as.vector(
            sapply(BinMethods, function(x){
                rep(x, length=nS)}))
    df <- data.frame(
        BinMethods=tmp,
        Samples=Label.Samples,
        FNR=out)
    df$BinMethods <- factor(df$BinMethods,
        levels=unique(df$BinMethods))
    df$Samples <- factor(df$Samples,
        levels=unique(df$Samples))
    df$FNR[which(df$FNR < 0)] <- 0
    df
}

# PR
Plot_PR <- function(){
    df <- aggregatePR()
    # Plot
    gg <- ggplot(df, aes(x=BinMethods, y=Samples, fill=PR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(Samples~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="PR", limits=c(0, max(df$PR)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=1, size=9))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave("plot/PR.png", gg, dpi=120, width=10, height=12)
}

aggregatePR <- function(){
    out <- c()
    nM <- length(BinMethods)
    nS <- length(Samples)
    for(i in 1:nM){
        for(j in 1:nS){
            inputfile <- paste0("output/",
                BinMethods[i], "/PR/", Samples[j], ".RData")
            load(inputfile)
            pr <- unlist(pr)
            target <- which(!is.nan(pr))
            out = c(out, mean(pr[target]))
        }
    }
    tmp <- as.vector(
            sapply(BinMethods, function(x){
                rep(x, length=nS)}))
    df <- data.frame(
        BinMethods=tmp,
        Samples=Label.Samples,
        PR=out)
    df$BinMethods <- factor(df$BinMethods,
        levels=unique(df$BinMethods))
    df$Samples <- factor(df$Samples,
        levels=unique(df$Samples))
    df$PR[which(df$PR < 0)] <- 0
    df
}

# TR
TR <- function(trueCaH){
    tc <- lapply(trueCaH, function(t){
        length(which(t == 1)) / length(t)
    })
    sum(unlist(tc))
}

Plot_TR <- function(){
    df <- aggregateTR()
    # Plot
    gg <- ggplot(df, aes(x=Samples, y=TR, fill=Samples))
    gg <- gg + geom_bar(stat="identity")
    gg <- gg + theme(axis.text.y=element_text(size = 30), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -90), axis.title.x=element_blank())
    gg <- gg + theme(legend.position="none")
    # Save
    ggsave("plot/TR.png", gg, dpi=240, width=6, height=9)
}

aggregateTR <- function(){
    out <- c()
    nS <- length(Samples)
    for(i in 1:nS){
        inputfile <- paste0("data/groundtruth/", Samples[i], ".RData")
        load(inputfile)
        out = c(out, TR(trueCaH))
    }
    df <- data.frame(
        Samples=Label.Samples,
        TR=out)
    df$Samples <- factor(df$Samples,
        levels=unique(df$Samples))
    df$TR[which(df$TR < 0)] <- 0
    df
}

# AUC
Plot_AUC <- function(){
    df <- aggregateAUC()
    # Plot
    gg <- ggplot(df, aes(x=Methods, y=Samples, fill=AUC))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(Samples~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="AUC", limits=c(0, 1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=1, size=9))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave("plot/AUC.png", gg, dpi=120, width=10, height=12)
}

aggregateAUC <- function(){
    out <- c()
    nM <- length(Methods)
    nS <- length(Samples)
    for(i in 1:nM){
        for(j in 1:nS){
            inputfile <- paste0("output/",
                Methods[i], "/AUC/", Samples[j], ".RData")
            load(inputfile)
            auc <- unlist(auc)
            target <- which(!is.nan(auc))
            out = c(out, mean(auc[target]))
        }
    }
    tmp <- as.vector(
            sapply(Methods, function(x){
                rep(x, length=nS)}))
    df <- data.frame(
        Methods=tmp,
        Samples=Label.Samples,
        AUC=out)
    df$Methods <- factor(df$Methods,
        levels=unique(df$Methods))
    df$Samples <- factor(df$Samples,
        levels=unique(df$Samples))
    df
}

# AUCPR
Plot_AUCPR <- function(){
    df <- aggregateAUCPR()
    # Plot
    gg <- ggplot(df, aes(x=Methods, y=Samples, fill=AUCPR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(Samples~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="AUCPR", limits=c(0, max(df$AUCPR)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=1, size=9))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave("plot/AUCPR.png", gg, dpi=120, width=10, height=12)
}

aggregateAUCPR <- function(){
    out <- c()
    nM <- length(Methods)
    nS <- length(Samples)
    for(i in 1:nM){
        for(j in 1:nS){
            inputfile <- paste0("output/",
                Methods[i], "/AUCPR/", Samples[j], ".RData")
            load(inputfile)
            aucpr <- unlist(aucpr)
            target <- which(!is.nan(aucpr))
            out = c(out, mean(aucpr[target]))
        }
    }
    tmp <- as.vector(
            sapply(Methods, function(x){
                rep(x, length=nS)}))
    df <- data.frame(
        Methods=tmp,
        Samples=Label.Samples,
        AUCPR=out)
    df$Methods <- factor(df$Methods,
        levels=unique(df$Methods))
    df$Samples <- factor(df$Samples,
        levels=unique(df$Samples))
    df
}

.ggdefault_cols <- function(n){
    hcl(h=seq(15, 375-360/n, length=n)%%360, c=100, l=65)
}

Plot_ROC_AUC_F <- function(roc, auc, f, outfile){
    colvec <- .ggdefault_cols(length(roc))
    png(file=outfile, width=1000, height=500)
    par(ps=22)
    par(plt=c(0.07, 0.93, 0.12, 0.88))
    for(i in seq_along(roc)){
        plot(roc[[i]],
            col=colvec[i],
            cex=2, pch=16, type="b", xlab="FPR", ylab="TPR",
            main="", xlim=c(0, 1), ylim=c(0, 1))
        par(new=TRUE)
    }
    # Legend
    l <- length(roc)
    # Human_Germline_Female, Mouse_Uterus, Human_FetalKidney
    if(1 <= l && l < 10){
        my.x <- 0.53
        my.y <- 0.62
        my.cex <- 1.5
    }
    # Human_NicotinehESCs_Nicotine
    if(11 <= l){
        my.x <- 0.8
        my.y <- 0.62
        my.cex <- 0.6
    }
    # Legend
    if(is.null(f)){
        legend.text <- paste0("AUC: ", unlist(auc))
    }else{
        legend.text <- paste0("AUC: ", unlist(auc),
            ", F: ", round(unlist(f), digits=7))
    }
    legend(my.x, my.y, legend.text,
        cex=my.cex, bty="n",
        col=colvec[seq_along(roc)], pch=16)
    dev.off()
}

Plot_PRC_AUCPR_MCC <- function(trueCaH, prc, aucpr, mcc, outfile){
    colvec <- .ggdefault_cols(length(prc))
    ymax <- max(unlist(lapply(prc, function(x){
        max(x$y, na.rm=TRUE)})))
    ylim=c(0, ymax)
    baseline <- lapply(trueCaH, function(x){
        length(which(x == 1)) / length(x)
    })
    png(file=outfile, width=1000, height=500)
    par(ps=22)
    par(plt=c(0.07, 0.93, 0.12, 0.88))
    for(i in seq_along(prc)){
        plot(prc[[i]],
            col=colvec[i],
            cex=2, pch=16, type="b", xlab="Recall (Sensitivity)",
            ylab="Precision (PPV)",
            main="", xlim=c(0, 1), ylim=ylim)
        abline(h = baseline[[i]], col=colvec[i])
        par(new=TRUE)
    }
    # Legend
    l <- length(prc)
    # Human_Germline_Female, Mouse_Uterus, Human_FetalKidney
    if(1 <= l && l < 10){
        my.x <- 0.53
        my.y <- 0.9 * ymax
        my.cex <- 1.5
    }
    # Human_NicotinehESCs_Nicotine
    if(11 <= l && l < 25){
        my.x <- 0.8
        my.y <- 0.9 * ymax
        my.cex <- 0.6
    }
    # Legend
    if(is.null(mcc)){
        legend.text <- paste0("AUCPR: ", unlist(aucpr))
    }else{
        legend.text <- paste0("AUCPR: ", unlist(aucpr),
            ", MCC: ", round(unlist(mcc), digits=7))
    }
    legend(my.x, my.y, legend.text,
        cex=my.cex, bty="n",
        col=colvec[seq_along(prc)], pch=16)
    dev.off()
}

multicciInfo <- function(lposition, rposition, ligandid, receptorid, out){
    # Setting
    colvec <- c(brewer.pal(9, "Set1"),
        brewer.pal(8, "Set2"),
        brewer.pal(8, "Dark2"))
    # L/R Pattern
    LPattern <- rep(0, length=dim(out$tnsr)[1])
    RPattern <- rep(0, length=dim(out$tnsr)[2])
    LPattern[lposition] <- 1
    RPattern[rposition] <- 1
    # CaH
    trueCaH <- array(0, dim=dim(out$tnsr))
    dimnames(trueCaH) <- list(
        unique(names(celltypes)),
        unique(names(celltypes)),
        out$pairname
        )
    target <- grep(paste(c(ligandid, receptorid), collapse="_")
            , out$pairname)
    if(length(target) != 0){
        trueCaH[lposition, rposition,
            grep(paste(c(ligandid, receptorid), collapse="_")
                , out$pairname)] <- 1
        trueCaH <- as.vector(trueCaH)
        return(list(
            LPattern=LPattern,
            RPattern=RPattern,
            trueCaH=trueCaH
            ))
    }else{
        stop("No L-R pair is found in the input tensor!!!")
    }
}

HVG <- function(data){
    data <- data + 1E-10
    # ED
    set.seed(123)
    ed <- t(t(data)/estimateSizeFactorsForMatrix(data))
    # Gene-wise mean
    m <- rowMeans(data)
    # Gene-wise variance
    v <- rowVars(as.matrix(data))
    # Gene-wise CV2
    cv2 <- v / m^2
    # Filtering of genes
    useForFit <- m >= unname(quantile(m[which( cv2 > 1)], .5 )) # & spikeins
    # Fitting
    set.seed(123)
    fit <- glmgam.fit(cbind(a0 = 1, a1tilde = 1/m[useForFit]), cv2[useForFit])
    # Parameters
    a0 <- unname(fit$coefficients["a0"])
    a1 <- unname(fit$coefficients["a1tilde"])
    afit <- a1 / m + a0
    varFitRatio <- v / (afit*m^2)
    varorder <- order(varFitRatio, decreasing=T)
    oed <- ed[varorder, ]
    xg <- exp(seq(min(log(m[m>0])), max(log(m)), length.out=1000))
    vfit <- a1/xg + a0
    df <- ncol(ed) - 1
    # p-valuse
    pval <- pchisq(varFitRatio*df, df=df, lower.tail=F)
    # Multiple test adjusting
    adj.pval <- p.adjust(pval, "fdr")
    # Output
    list(ed=ed, m=m, v=v, cv2=cv2, useForFit=useForFit, fit=fit, a0=a0, a1=a1, afit=afit, varFitRatio=varFitRatio, varorder=varorder, oed=oed, xg=xg, vfit=vfit, df=df, pval=pval, adj.pval=adj.pval)
}

ggdefault_cols <- function(n){
    hcl(h=seq(15, 375-360/n, length=n)%%360, c=100, l=65)
}

.lapply_pb <- function(X, FUN, ...)
{
 env <- environment()
 pb_Total <- length(X)
 counter <- 0
 pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

 # wrapper around FUN
 wrapper <- function(...){
   curVal <- get("counter", envir = env)
   assign("counter", curVal +1 ,envir=env)
   setTxtProgressBar(get("pb", envir=env), curVal +1)
   FUN(...)
 }
 res <- lapply(X, wrapper, ...)
 close(pb)
 res
}

AUC <- function(score, trueCaH){
    pred <- prediction(score, trueCaH)
    round(performance(pred, "auc")@y.values[[1]], 2)
}

AUCPR <- function(score, trueCaH){
    pred <- prediction(score, trueCaH)
    round(performance(pred, "aucpr")@y.values[[1]], 2)
}

Fmeasure <- function(predict, label){
    TP <- length(intersect(which(predict == 1), which(label == 1)))
    FP <- length(intersect(which(predict == 1), which(label == 0)))
    FN <- length(intersect(which(predict == 0), which(label == 1)))
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    (2 * Recall * Precision) / (Recall + Precision)
}

MCC <- function(predict, label){
    cor(predict, label, method="pearson")
}

FPR <- function(predict, label){
    FP <- length(intersect(which(predict == 1), which(label == 0)))
    TN <- length(intersect(which(predict == 0), which(label == 0)))
    FP / (TN + FP)
}

FNR <- function(predict, label){
    TP <- length(intersect(which(predict == 1), which(label == 1)))
    FN <- length(intersect(which(predict == 0), which(label == 1)))
    FN / (TP + FN)
}

PR <- function(predict){
    length(which(predict == 1)) / length(predict)
}

Tensor2Vec <- function(ncelltypes, cci, out){
    counter <- 1
    score <- rep(0, length=ncelltypes^2*cci$nPair)
    for(j in seq_len(cci$nPair)){
        score[counter:(counter+ncelltypes^2-1)] <- as.vector(out[,,j])
        counter <- counter + ncelltypes^2
    }
    score
}

ROCCurve = function(score, actual){
    pred <- prediction(score, actual)
    nF <- sum(actual == 0)
    nT <- sum(actual == 1)
    x <- pred@fp[[1]] / nF
    y <- pred@tp[[1]] / nT
    list(x=x, y=y)
}

PRCCurve <- function(score, actual){
    pred <- prediction(score, actual)
    perf <- performance(pred, "prec", "rec")
    x <- perf@x.values[[1]]
    y <- perf@y.values[[1]]
    list(x=x, y=y)
}

ROC_AUC_BIN_F <- function(ncelltypes, cif, trueCaH, out, outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9, outfile10, pval=FALSE){
    # Scoring
    if(pval){
        score <- Tensor2Vec(ncelltypes, cif, 1 - out$pval)
        bin <- 1 - score
	    target <- which(p.adjust(bin) < 0.1)
	    # Binarization
	    bin[target] <- 1
	    bin[setdiff(seq(length(bin)), target)] <- 0
	    # F-measure
        f <- lapply(trueCaH, function(t, bin){
		    Fmeasure(bin, t)
    	}, bin=bin)
        # MCC
        mcc <- lapply(trueCaH, function(t, bin){
            MCC(bin, t)
        }, bin=bin)
        # FPR
        fpr <- lapply(trueCaH, function(t, bin){
            FPR(bin, t)
        }, bin=bin)
        # FNR
        fnr <- lapply(trueCaH, function(t, bin){
            FNR(bin, t)
        }, bin=bin)
        # PR
        pr <- PR(bin)
    }else{
        score <- Tensor2Vec(ncelltypes, cif, out$tnsr)
        bin <- NULL
        f <- NULL
        mcc <- NULL
        fpr <- NULL
        fnr <- NULL
        pr <- NULL
    }
    # AUC
    auc <- lapply(trueCaH, function(t, score){
        AUC(score, t)
    }, score=score)
    # ROC
    roc <- lapply(trueCaH, function(t, score){
        ROCCurve(score, t)
    }, score=score)
    # AUCPR
    aucpr <- lapply(trueCaH, function(t, score){
        AUCPR(score, t)
    }, score=score)
    # PRC
    prc <- lapply(trueCaH, function(t, score){
        PRCCurve(score, t)
    }, score=score)
    # Save
    save(roc, file=outfile1)
    save(auc, file=outfile2)
    save(bin, file=outfile3)
    save(f, file=outfile4)
    save(prc, file=outfile5)
    save(aucpr, file=outfile6)
    save(mcc, file=outfile7)
    save(fpr, file=outfile8)
    save(fnr, file=outfile9)
    save(pr, file=outfile10)
}

BIN <- function(out, x){
	# L <- scTensor:::.HCLUST(out$ligand[x[1],])
	# R <- scTensor:::.HCLUST(out$receptor[x[2],])
	# LR <- scTensor:::.HCLUST(out$lrpair[x[3],])
    # L <- .HCLUST(out$ligand[x[1],])
    # R <- .HCLUST(out$receptor[x[2],])
    # LR <- .HCLUST(out$lrpair[x[3],])
    # L <- .DBSCAN(out$ligand[x[1],])
    # R <- .DBSCAN(out$receptor[x[2],])
    # LR <- .DBSCAN(out$lrpair[x[3],])
    # L <- .KMEANS(out$ligand[x[1],])
    # R <- .KMEANS(out$receptor[x[2],])
    # LR <- .KMEANS(out$lrpair[x[3],])
    # L <- .GMM(out$ligand[x[1],])
    # R <- .GMM(out$receptor[x[2],])
    # LR <- .GMM(out$lrpair[x[3],])
    # L <- .SD(out$ligand[x[1],])
    # R <- .SD(out$receptor[x[2],])
    # LR <- .SD(out$lrpair[x[3],])
    L <- .MAD(out$ligand[x[1],])
    R <- .MAD(out$receptor[x[2],])
    LR <- .MAD(out$lrpair[x[3],])
	L[which(L == "selected")] <- 1
	R[which(R == "selected")] <- 1
	LR[which(LR == "selected")] <- 1
	L[which(L == "not selected")] <- 0
	R[which(R == "not selected")] <- 0
	LR[which(LR == "not selected")] <- 0
	L <- as.numeric(L)
	R <- as.numeric(R)
	LR <- as.numeric(LR)
	list(ligand=L, receptor=R, lrpair=LR)
}

ROCAUCF <- function(trueCaH, score, bin){
	lapply(trueCaH, function(t, other){
	    	score <- other$score
	    	bin <- other$bin
	        tmp <- unlist(lapply(score, function(s, t){
	            AUC(as.vector(s@data), t)
	        }, t=t))
	        maxposition <- which(max(tmp) == tmp)[1]
	        # AUC
	        auc <- max(tmp)
	        # ROC
            roc <- ROCCurve(as.vector(score[[maxposition]]@data), t)
	        # F-measure
	        f <- Fmeasure(as.vector(bin[[maxposition]]@data), t)
	        list(auc=auc, roc=roc, f=f)
	    }, other=list(score=score, bin=bin))
}

ROCAUCF2 <- function(trueCaH, score, bin){
    lapply(trueCaH, function(t, other){
            # Max ROC-AUC Search
            score <- other$score
            bin <- other$bin
            tmp <- unlist(lapply(score, function(s, t){
                AUC(as.vector(s@data), t)
            }, t=t))
            maxposition <- which(max(tmp) == tmp)[1]
            # AUC
            auc <- max(tmp)
            # ROC
            roc <- ROCCurve(as.vector(score[[maxposition]]@data), t)
            # F-measure
            f <- Fmeasure(as.vector(bin[[maxposition]]@data), t)
            # Max ROC-AUC Search 2
            tmp2 <- unlist(lapply(score, function(s, t){
                AUCPR(as.vector(s@data), t)
            }, t=t))
            maxposition2 <- which(max(tmp2) == tmp2)[1]
            # AUC
            aucpr <- max(tmp2)
            # PRC
            prc <- PRCCurve(as.vector(score[[maxposition2]]@data), t)
            # MCC
            mcc <- MCC(as.vector(bin[[maxposition2]]@data), t)
            # FPR
            fpr <- FPR(as.vector(bin[[maxposition2]]@data), t)
            # FNR
            fnr <- FNR(as.vector(bin[[maxposition2]]@data), t)
            # PR
            pr <- PR(as.vector(bin[[maxposition2]]@data))
            # Output
            list(auc=auc, roc=roc, f=f, aucpr=aucpr, prc=prc,
                mcc=mcc, fpr=fpr, fnr=fnr, pr=pr)
        }, other=list(score=score, bin=bin))
}

ROC_AUC_BIN_F_previous_sctensor <- function(trueCaH, out, outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9, outfile10){
    # Scoring
    score <- apply(out$index, 1, function(x){
        nnTensor::recTensor(x[4],
            list(as.matrix(out$ligand[x[1], ]),
                as.matrix(out$receptor[x[2], ]),
                as.matrix(out$lrpair[x[3], ])),
            		reverse=TRUE)
    })
    # Binarization
    bin <- apply(out$index, 1, function(x){
    	tmp <- BIN(out, x)
    	ligand <- tmp$ligand
    	receptor <- tmp$receptor
    	lrpair <- tmp$lrpair
        nnTensor::recTensor(1,
            list(as.matrix(ligand),
                as.matrix(receptor),
                as.matrix(lrpair)), reverse=TRUE)
    })
    # ROC / AUC
    rocaucf <- ROCAUCF2(trueCaH, score, bin)
    roc <- lapply(rocaucf, function(raf){
    	raf$roc
    })
    auc <- lapply(rocaucf, function(raf){
    	raf$auc
    })
    f <- lapply(rocaucf, function(raf){
    	raf$f
    })
    prc <- lapply(rocaucf, function(raf){
        raf$prc
    })
    aucpr <- lapply(rocaucf, function(raf){
        raf$aucpr
    })
    mcc <- lapply(rocaucf, function(raf){
        raf$mcc
    })
    fpr <- lapply(rocaucf, function(raf){
        raf$fpr
    })
    fnr <- lapply(rocaucf, function(raf){
        raf$fnr
    })
    pr <- lapply(rocaucf, function(raf){
        raf$pr
    })
    # Save
    save(roc, file=outfile1)
    save(auc, file=outfile2)
    save(bin, file=outfile3)
    save(f, file=outfile4)
    save(prc, file=outfile5)
    save(aucpr, file=outfile6)
    save(mcc, file=outfile7)
    save(fpr, file=outfile8)
    save(fnr, file=outfile9)
    save(pr, file=outfile10)
}

# previous_sctensorはよくなったが、sctensorはよくならなかった
.DBSCAN <- function(x){
    cluster <- .searchEPS(x)
    max1 <- max(x[which(cluster == 1)])
    max2 <- max(x[which(cluster == 2)])
    if(max1 > max2){
        cluster[which(cluster == 1)] <- "selected"
        cluster[which(cluster == 2)] <- "not selected"
    }else{
        cluster[which(cluster == 1)] <- "not selected"
        cluster[which(cluster == 2)] <- "selected"
    }
    cluster
}

# single, complete, ward.D2, average, median, centroid, mcquittyはやった
# medianが良さそう
.HCLUST <- function(x){
    out <- hclust(dist(x), method="mcquitty")
    cluster <- cutree(out, 2)
    max1 <- max(x[which(cluster == 1)])
    max2 <- max(x[which(cluster == 2)])
    if(max1 > max2){
        cluster[which(cluster == 1)] <- "selected"
        cluster[which(cluster == 2)] <- "not selected"
    }else{
        cluster[which(cluster == 1)] <- "not selected"
        cluster[which(cluster == 2)] <- "selected"
    }
    cluster
}

.KMEANS <- function(x){
    cluster <- kmeans(x, centers=2)$cluster
    max1 <- max(x[which(cluster == 1)])
    max2 <- max(x[which(cluster == 2)])
    if(max1 > max2){
        cluster[which(cluster == 1)] <- "selected"
        cluster[which(cluster == 2)] <- "not selected"
    }else{
        cluster[which(cluster == 1)] <- "not selected"
        cluster[which(cluster == 2)] <- "selected"
    }
    cluster
}

.GMM <- function(x){
    cluster <- Mclust(x, 2)$classification
    max1 <- max(x[which(cluster == 1)])
    max2 <- max(x[which(cluster == 2)])
    if(max1 > max2){
        cluster[which(cluster == 1)] <- "selected"
        cluster[which(cluster == 2)] <- "not selected"
    }else{
        cluster[which(cluster == 1)] <- "not selected"
        cluster[which(cluster == 2)] <- "selected"
    }
    cluster
}

.SD <- function(x, thr=1){
    cluster <- x - mean(x) / sd(x)
    target <- which(cluster >= thr)
    cluster[target] <- "selected"
    cluster[setdiff(seq(cluster), target)] <- "not selected"
    cluster
}

.MAD <- function(x, thr=1.0){
    cluster <- abs(x - median(x))
    target <- which((x - median(x)) >= thr*median(cluster))
    cluster[target] <- "selected"
    cluster[setdiff(seq(cluster), target)] <- "not selected"
    cluster
}

.searchEPS <- function(x){
    div <- 100:1
    for(i in seq(div)){
        cls <- dbscan(as.matrix(x),
            eps=(max(x)-min(x))/div[i],
            minPts=length(x)/10)$cluster
        if(length(unique(cls)) == 2){
            break
        }
    }
    if(length(unique(cls)) != 2){
        tc <- table(cls)
        maxpos <- which(tc == max(tc))[1]
        target <- which(cls == maxpos)
        cls[target] <- 1
        cls[setdiff(seq(cls), target)] <- 2
    }
    cls
}

BIN_2 <- function(out, x){
    # L <- scTensor:::.HCLUST(out$ligand[x[1],])
    # R <- scTensor:::.HCLUST(out$receptor[x[2],])
    # L <- .HCLUST(out$ligand[x[1],])
    # R <- .HCLUST(out$receptor[x[2],])
    # L <- .DBSCAN(out$ligand[x[1],])
    # R <- .DBSCAN(out$receptor[x[2],])
    # L <- .KMEANS(out$ligand[x[1],])
    # R <- .KMEANS(out$receptor[x[2],])
    # L <- .GMM(out$ligand[x[1],])
    # R <- .GMM(out$receptor[x[2],])
    # L <- .SD(out$ligand[x[1],])
    # R <- .SD(out$receptor[x[2],])
    L <- .MAD(out$ligand[x[1],])
    R <- .MAD(out$receptor[x[2],])
    L[which(L == "selected")] <- 1
    R[which(R == "selected")] <- 1
    L[which(L == "not selected")] <- 0
    R[which(R == "not selected")] <- 0
    L <- as.numeric(L)
    R <- as.numeric(R)
    list(ligand=L, receptor=R)
}

ROC_AUC_BIN_F_sctensor <- function(trueCaH, out, outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9, outfile10){
    # Scoring
    score <- apply(out$index, 1, function(x){
        nnTensor::recTensor(x[4],
            list(as.matrix(out$ligand[x[1], ]),
                as.matrix(out$receptor[x[2], ]),
                as.matrix(out$lrpair[x[1], x[2], ]@data)), reverse=TRUE)
    })
    # Binarization
    bin.core <- out$lrpair@data
    for(i in seq(dim(bin.core)[3])){
        # tmp <- scTensor:::.HCLUST(as.vector(bin.core[,,i]))
        # tmp <- .HCLUST(as.vector(bin.core[,,i]))
        # tmp <- .DBSCAN(as.vector(bin.core[,,i]))
        # tmp <- .KMEANS(as.vector(bin.core[,,i]))
        # tmp <- .GMM(as.vector(bin.core[,,i]))
        # tmp <- .SD(as.vector(bin.core[,,i]))
        tmp <- .MAD(as.vector(bin.core[,,i]))
        tmp[which(tmp == "selected")] <- 1
        tmp[which(tmp == "not selected")] <- 0
        tmp <- as.numeric(tmp)
        dim(tmp) <- dim(bin.core)[1:2]
        bin.core[,,i] <- tmp

    }
    bin <- apply(out$index, 1, function(x){
        tmp <- BIN_2(out, x)
        ligand <- tmp$ligand
        receptor <- tmp$receptor
        lrpair <- bin.core[x[1],x[2],]
        nnTensor::recTensor(1,
            list(as.matrix(ligand),
                as.matrix(receptor),
                as.matrix(lrpair)), reverse=TRUE)
    })
    # ROC / AUC
    rocaucf <- ROCAUCF2(trueCaH, score, bin)
    roc <- lapply(rocaucf, function(raf){
    	raf$roc
    })
    auc <- lapply(rocaucf, function(raf){
    	raf$auc
    })
    f <- lapply(rocaucf, function(raf){
    	raf$f
    })
    prc <- lapply(rocaucf, function(raf){
        raf$prc
    })
    aucpr <- lapply(rocaucf, function(raf){
        raf$aucpr
    })
    mcc <- lapply(rocaucf, function(raf){
        raf$mcc
    })
    fpr <- lapply(rocaucf, function(raf){
        raf$fpr
    })
    fnr <- lapply(rocaucf, function(raf){
        raf$fnr
    })
    pr <- lapply(rocaucf, function(raf){
        raf$pr
    })
    # Save
    save(roc, file=outfile1)
    save(auc, file=outfile2)
    save(bin, file=outfile3)
    save(f, file=outfile4)
    save(prc, file=outfile5)
    save(aucpr, file=outfile6)
    save(mcc, file=outfile7)
    save(fpr, file=outfile8)
    save(fnr, file=outfile9)
    save(pr, file=outfile10)
}

groundTruth <- function(spl){
	cif <- cciInfo[[spl]]
	target <- length(grep("CCI", names(cif)))
	ncelltypes <- length(cif$CCI1$LPattern)
	trueCaH <- list()
	length(trueCaH) <- target
	counter <- 1
	# Each CCI1, 2, 3, ...
	for(j in seq_along(target)){
	    L <- eval(parse(text=paste0("cif$CCI", j, "$LPattern")))
	    R <- eval(parse(text=paste0("cif$CCI", j, "$RPattern")))
	    nGene <- eval(parse(text=paste0("cif$CCI", j, "$nGene")))
	    # trueCaH
	    tmpTrueCaH <- rep(as.vector(base::outer(L, R, "*")), nGene)
	    # Setting true CaH as 1
	    trueCaH[[j]] <- rep(0, length=ncelltypes^2*cif$nPair)
	    trueCaH[[j]][counter:(counter+length(tmpTrueCaH)-1)] <- tmpTrueCaH

	    counter <- counter + length(tmpTrueCaH)
	}
	list(trueCaH=trueCaH, ncelltypes=ncelltypes, cif=cif, spl=spl)
}

cciInfo <- list(
    ###### 1. 3 Celltypes, 1 CCI Patterns, One-by-one
    "3Celltypes_1CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0),
         RPattern=c(0,1,0),
         nGene=50,
         fc="E2")
      ),
    ###### 2. 5 Celltypes, 1 CCI Patterns, One-by-one
    "5Celltypes_1CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0),
         RPattern=c(0,1,0,0,0),
         nGene=50,
         fc="E2")
      ),
    ###### 3. 10 Celltypes, 1 CCI Patterns, One-by-one
    "10Celltypes_1CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,1,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2")
      ),
    ###### 4. 20 Celltypes, 1 CCI Patterns, One-by-one
    "20Celltypes_1CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2")
      ),
    ###### 5. 30 Celltypes, 1 CCI Patterns, One-by-one
    "30Celltypes_1CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2")
      ),

    ###### 6. 3 Celltypes, 3 CCI Patterns, One-by-one
    # Fig.のCase I
    "3Celltypes_3CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0),
         RPattern=c(0,1,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,1,0),
         RPattern=c(0,0,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(0,0,1),
         RPattern=c(1,0,0),
         nGene=50,
         fc="E2")
      ),
    ###### 7. 5 Celltypes, 3 CCI Patterns, One-by-one
    "5Celltypes_3CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0),
         RPattern=c(0,1,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,1,0,0,0),
         RPattern=c(0,0,1,0,0),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(0,0,1,0,0),
         RPattern=c(0,0,0,1,0),
         nGene=50,
         fc="E2")
      ),
    ###### 8. 10 Celltypes, 3 CCI Patterns, One-by-one
    "10Celltypes_3CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,1,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,1,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,1,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(0,0,1,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,1,0,0,0,0,0,0),
         nGene=50,
         fc="E2")
      ),
    ###### 9. 20 Celltypes, 3 CCI Patterns, One-by-one
    "20Celltypes_3CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2")
      ),
    ###### 10. 30 Celltypes, 3 CCI Patterns, One-by-one
    "30Celltypes_3CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2")
      ),

    ###### 11. 3 Celltypes, 5 CCI Patterns, One-by-one
    "3Celltypes_5CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0),
         RPattern=c(0,1,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,1,0),
         RPattern=c(0,0,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(0,0,1),
         RPattern=c(1,0,0),
         nGene=50,
         fc="E2"),
      CCI4=list(
         LPattern=c(0,1,0),
         RPattern=c(1,0,0),
         nGene=50,
         fc="E2"),
      CCI5=list(
         LPattern=c(0,0,1),
         RPattern=c(0,1,0),
         nGene=50,
         fc="E2")
      ),
    ###### 12. 5 Celltypes, 5 CCI Patterns, One-by-one
    "5Celltypes_5CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0),
         RPattern=c(0,1,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,1,0,0,0),
         RPattern=c(0,0,1,0,0),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(0,0,1,0,0),
         RPattern=c(0,0,0,1,0),
         nGene=50,
         fc="E2"),
      CCI4=list(
         LPattern=c(0,0,0,1,0),
         RPattern=c(0,0,0,0,1),
         nGene=50,
         fc="E2"),
      CCI5=list(
         LPattern=c(0,0,0,0,1),
         RPattern=c(1,0,0,0,0),
         nGene=50,
         fc="E2")
      ),
    ###### 13. 10 Celltypes, 5 CCI Patterns, One-by-one
    "10Celltypes_5CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,1,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,1,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,1,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(0,0,1,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,1,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI4=list(
         LPattern=c(0,0,0,1,0,0,0,0,0,0),
         RPattern=c(0,0,0,0,1,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI5=list(
         LPattern=c(0,0,0,0,1,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,0,0,0,0),
         nGene=50,
         fc="E2")
      ),
    ###### 14. 20 Celltypes, 5 CCI Patterns, One-by-one
    "20Celltypes_5CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI4=list(
         LPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI5=list(
         LPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2")
      ),
    ###### 15. 30 Celltypes, 5 CCI Patterns, One-by-one
    "30Celltypes_5CCIPatterns_OnetoOne"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI4=list(
         LPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI5=list(
         LPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2")
      ),

    ###### 16. 3 Celltypes, 1 CCI Patterns, Many-by-many
    "3Celltypes_1CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0),
         RPattern=c(0,1,1),
         nGene=50,
         fc="E2")
      ),
    ###### 17. 5 Celltypes, 1 CCI Patterns, Many-by-many
    "5Celltypes_1CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,0,0,0),
         RPattern=c(0,0,1,1,1),
         nGene=50,
         fc="E2")
      ),
    ###### 18. 10 Celltypes, 1 CCI Patterns, Many-by-many
    "10Celltypes_1CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,1,1,1,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,1,1,1,1),
         nGene=50,
         fc="E2")
      ),
    ###### 19. 20 Celltypes, 1 CCI Patterns, Many-by-many
    "20Celltypes_1CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2")
      ),
    ###### 20. 30 Celltypes, 1 CCI Patterns, Many-by-many
    "30Celltypes_1CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2")
      ),

    ###### 21. 3 Celltypes, 3 CCI Patterns, Many-by-many
    # Fig.のCase II
    "3Celltypes_3CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0),
         RPattern=c(0,1,1),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,0,1),
         RPattern=c(1,1,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(1,1,1),
         RPattern=c(1,0,1),
         nGene=50,
         fc="E2")
      ),
    ###### 22. 5 Celltypes, 3 CCI Patterns, Many-by-many
    "5Celltypes_3CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,0,0,0),
         RPattern=c(0,0,1,1,1),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,0,0,0,1),
         RPattern=c(1,1,1,1,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(1,1,1,1,1),
         RPattern=c(1,1,0,0,1),
         nGene=50,
         fc="E2")
      ),
    ###### 23. 10 Celltypes, 3 CCI Patterns, Many-by-many
    "10Celltypes_3CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,1,1,1,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,1,1,1,1),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,0,0,0,0,0,0,0,0,1),
         RPattern=c(1,1,1,1,1,1,1,1,1,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(1,1,1,1,1,1,1,1,1,1),
         RPattern=c(1,1,0,0,0,0,0,0,1,1),
         nGene=50,
         fc="E2")
      ),
    ###### 24. 20 Celltypes, 3 CCI Patterns, Many-by-many
    "20Celltypes_3CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
         RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
         RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
         nGene=50,
         fc="E2")
      ),
    ###### 25. 30 Celltypes, 3 CCI Patterns, Many-by-many
    "30Celltypes_3CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
         RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
         RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
         nGene=50,
         fc="E2")
      ),

    ###### 26. 3 Celltypes, 5 CCI Patterns, Many-by-many
    "3Celltypes_5CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,0,0),
         RPattern=c(0,1,1),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,0,1),
         RPattern=c(1,1,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(1,1,1),
         RPattern=c(1,0,1),
         nGene=50,
         fc="E2"),
      CCI4=list(
         LPattern=c(1,0,1),
         RPattern=c(0,1,0),
         nGene=50,
         fc="E2"),
      CCI5=list(
         LPattern=c(0,1,0),
         RPattern=c(1,0,1),
         nGene=50,
         fc="E2")
      ),
    ###### 27. 5 Celltypes, 5 CCI Patterns, Many-by-many
    "5Celltypes_5CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,0,0,0),
         RPattern=c(0,0,1,1,1),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,0,0,0,1),
         RPattern=c(1,1,1,1,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(1,1,1,1,1),
         RPattern=c(1,1,0,0,1),
         nGene=50,
         fc="E2"),
      CCI4=list(
         LPattern=c(1,0,1,0,1),
         RPattern=c(0,1,0,1,0),
         nGene=50,
         fc="E2"),
      CCI5=list(
         LPattern=c(0,1,0,1,0),
         RPattern=c(1,0,1,0,1),
         nGene=50,
         fc="E2")
      ),
    ###### 28. 10 Celltypes, 5 CCI Patterns, Many-by-many
    "10Celltypes_5CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,1,1,1,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,1,1,1,1),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,0,0,0,0,0,0,0,0,1),
         RPattern=c(1,1,1,1,1,1,1,1,1,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(1,1,1,1,1,1,1,1,1,1),
         RPattern=c(1,1,0,0,0,0,0,0,1,1),
         nGene=50,
         fc="E2"),
      CCI4=list(
         LPattern=c(1,0,1,0,1,0,1,0,1,0),
         RPattern=c(0,1,0,1,0,1,0,1,0,1),
         nGene=50,
         fc="E2"),
      CCI5=list(
         LPattern=c(0,1,0,1,0,1,0,1,0,1),
         RPattern=c(1,0,1,0,1,0,1,0,1,0),
         nGene=50,
         fc="E2")
      ),
    ###### 29. 20 Celltypes, 5 CCI Patterns, Many-by-many
    "20Celltypes_5CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
         RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
         RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
         nGene=50,
         fc="E2"),
      CCI4=list(
         LPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
         RPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
         nGene=50,
         fc="E2"),
      CCI5=list(
         LPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
         RPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
         nGene=50,
         fc="E2")
      ),
    ###### 30. 30 Celltypes, 5 CCI Patterns, Many-by-many
    "30Celltypes_5CCIPatterns_ManytoMany"=
    list(nPair=1000,
      CCI1=list(
         LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
         nGene=50,
         fc="E2"),
      CCI2=list(
         LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
         RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
         nGene=50,
         fc="E2"),
      CCI3=list(
         LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
         RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
         nGene=50,
         fc="E2"),
      CCI4=list(
         LPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
         RPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
         nGene=50,
         fc="E2"),
      CCI5=list(
         LPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
         RPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
         nGene=50,
         fc="E2")
      )
)

nCells <- list(
    "3Celltypes_1CCIPatterns_OnetoOne"=rep(50,3),
    "5Celltypes_1CCIPatterns_OnetoOne"=rep(50,5),
    "10Celltypes_1CCIPatterns_OnetoOne"=rep(50,10),
    "20Celltypes_1CCIPatterns_OnetoOne"=rep(50,20),
    "30Celltypes_1CCIPatterns_OnetoOne"=rep(50,30),
    "3Celltypes_3CCIPatterns_OnetoOne"=rep(50,3),
    "5Celltypes_3CCIPatterns_OnetoOne"=rep(50,5),
    "10Celltypes_3CCIPatterns_OnetoOne"=rep(50,10),
    "20Celltypes_3CCIPatterns_OnetoOne"=rep(50,20),
    "30Celltypes_3CCIPatterns_OnetoOne"=rep(50,30),
    "3Celltypes_5CCIPatterns_OnetoOne"=rep(50,3),
    "5Celltypes_5CCIPatterns_OnetoOne"=rep(50,5),
    "10Celltypes_5CCIPatterns_OnetoOne"=rep(50,10),
    "20Celltypes_5CCIPatterns_OnetoOne"=rep(50,20),
    "30Celltypes_5CCIPatterns_OnetoOne"=rep(50,30),
    "3Celltypes_1CCIPatterns_ManytoMany"=rep(50,3),
    "5Celltypes_1CCIPatterns_ManytoMany"=rep(50,5),
    "10Celltypes_1CCIPatterns_ManytoMany"=rep(50,10),
    "20Celltypes_1CCIPatterns_ManytoMany"=rep(50,20),
    "30Celltypes_1CCIPatterns_ManytoMany"=rep(50,30),
    "3Celltypes_3CCIPatterns_ManytoMany"=rep(50,3),
    "5Celltypes_3CCIPatterns_ManytoMany"=rep(50,5),
    "10Celltypes_3CCIPatterns_ManytoMany"=rep(50,10),
    "20Celltypes_3CCIPatterns_ManytoMany"=rep(50,20),
    "30Celltypes_3CCIPatterns_ManytoMany"=rep(50,30),
    "3Celltypes_5CCIPatterns_ManytoMany"=rep(50,3),
    "5Celltypes_5CCIPatterns_ManytoMany"=rep(50,5),
    "10Celltypes_5CCIPatterns_ManytoMany"=rep(50,10),
    "20Celltypes_5CCIPatterns_ManytoMany"=rep(50,20),
    "30Celltypes_5CCIPatterns_ManytoMany"=rep(50, 30)
)
