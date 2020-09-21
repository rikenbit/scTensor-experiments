library("scTensor")
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
library("LRBase.Hsa.eg.db")
library("LRBase.Mmu.eg.db")
library("ggplot2")
library("viridis")
library("dbscan")

# Check 2.0.0
if(package.version("LRBase.Hsa.eg.db") != "2.0.0"){
    stop("LRBase.Hsa.eg.db is not 2.0.0!")
}
if(package.version("LRBase.Mmu.eg.db") != "2.0.0"){
    stop("LRBase.Mmu.eg.db is not 2.0.0!")
}

# For ROC, AUC, Binarization, F-measure
tensor.methods <- c(
    "labelpermutation_tensor", "labelpermutation2_tensor",
    "halpern_tensor", "cabelloaguilar_tensor")
pval.methods <- c(
    "labelpermutation", "labelpermutation2",
    "halpern", "cabelloaguilar")

# Time
Plot_Memory <- function(){
    df <- aggregateMemory()
    # Plot
    gg <- ggplot(df, aes(x=Methods, y=Samples, fill=GB))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(Samples~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="GB", limits=c(0, max(df$GB)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=13))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    outfile <- "plot/Memory.png"
    ggsave(outfile, gg, dpi=120, width=10, height=12)
}

aggregateMemory <- function(){
    out <- c()
    Methods = c("labelpermutation_tensor",
    "labelpermutation",
    "labelpermutation2_tensor",
    "labelpermutation2",
    "halpern_tensor",
    "halpern",
    "cabelloaguilar_tensor",
    "cabelloaguilar",
    "previous_sctensor",
    "sctensor")
    Samples = c("Human_FetalKidney",
        "Human_NicotinehESCs_Nicotine",
        "Human_Germline_Female",
        "Mouse_Uterus")
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
        Samples=c("FetalKidney", "NicotinehESCs", "Germline", "Uterus"),
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
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=13))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    outfile <- "plot/Time.png"
    ggsave(outfile, gg, dpi=120, width=10, height=12)
}

aggregateTime <- function(){
    out <- c()
    Methods = c("labelpermutation_tensor",
    "labelpermutation",
    "labelpermutation2_tensor",
    "labelpermutation2",
    "halpern_tensor",
    "halpern",
    "cabelloaguilar_tensor",
    "cabelloaguilar",
    "previous_sctensor",
    "sctensor")
    Samples = c("Human_FetalKidney",
        "Human_NicotinehESCs_Nicotine",
        "Human_Germline_Female",
        "Mouse_Uterus")
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
        Samples=c("FetalKidney", "NicotinehESCs", "Germline", "Uterus"),
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
    gg <- ggplot(df, aes(x=Methods, y=Samples, fill=F))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(Samples~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="F-measure", limits=c(0, max(df$F)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=13))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    outfile <- "plot/F.png"
    ggsave(outfile, gg, dpi=120, width=10, height=12)
}

aggregateF <- function(){
    out <- c()
    Methods = c("labelpermutation_tensor",
    "labelpermutation",
    "labelpermutation2_tensor",
    "labelpermutation2",
    "halpern_tensor",
    "halpern",
    "cabelloaguilar_tensor",
    "cabelloaguilar",
    "previous_sctensor",
    "sctensor")
    Samples = c("Human_FetalKidney",
        "Human_NicotinehESCs_Nicotine",
        "Human_Germline_Female",
        "Mouse_Uterus")
    nM <- length(Methods)
    nS <- length(Samples)
    for(i in 1:nM){
        for(j in 1:nS){
            inputfile <- paste0("output/",
                Methods[i], "/F/", Samples[j], ".RData")
            load(inputfile)
            f <- unlist(f)
            target <- which(!is.nan(f))
            out = c(out, mean(f[target]))
        }
    }
    tmp <- as.vector(
            sapply(Methods, function(x){
                rep(x, length=nS)}))
    df <- data.frame(
        Methods=tmp,
        Samples=c("FetalKidney", "NicotinehESCs", "Germline", "Uterus"),
        F=out)
    df$Methods <- factor(df$Methods,
        levels=unique(df$Methods))
    df$Samples <- factor(df$Samples,
        levels=unique(df$Samples))
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
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=13))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0, angle = -60), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    outfile <- "plot/AUC.png"
    ggsave(outfile, gg, dpi=120, width=10, height=12)
}

aggregateAUC <- function(){
    out <- c()
    Methods = c("labelpermutation_tensor",
    "labelpermutation",
    "labelpermutation2_tensor",
    "labelpermutation2",
    "halpern_tensor",
    "halpern",
    "cabelloaguilar_tensor",
    "cabelloaguilar",
    "previous_sctensor",
    "sctensor")
    Samples = c("Human_FetalKidney",
        "Human_NicotinehESCs_Nicotine",
        "Human_Germline_Female",
        "Mouse_Uterus")
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
        Samples=c("FetalKidney", "NicotinehESCs", "Germline", "Uterus"),
        AUC=out)
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

Fmeasure <- function(predict, label){
    TP <- length(intersect(which(predict == 1), which(label == 1)))
    FP <- length(intersect(which(predict == 1), which(label == 0)))
    FN <- length(intersect(which(predict == 0), which(label == 1)))
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    (2 * Recall * Precision) / (Recall + Precision)
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

ROC_AUC_BIN_F <- function(ncelltypes, cif, trueCaH, out, outfile1, outfile2, outfile3, outfile4, pval=FALSE){
    # Scoring
    if(pval){
        score <- Tensor2Vec(ncelltypes, cif, 1 - out$pval)
        bin <- score
	    target <- which(bin < 0.05)
	    # Binarization
	    bin[target] <- 0
	    bin[setdiff(seq(length(bin)), target)] <- 1
	    # F-measure
        f <- lapply(trueCaH, function(t, bin){
		    Fmeasure(bin, t)
    	}, bin=bin)
    }else{
        score <- Tensor2Vec(ncelltypes, cif, out$tnsr)
        bin <- NULL
        f <- NULL
    }
    # AUC
    auc <- lapply(trueCaH, function(t, score){
        AUC(score, t)
    }, score=score)
    # ROC
    roc <- lapply(trueCaH, function(t, score){
        ROCCurve(score, t)
    }, score=score)
    # Save
    save(roc, file=outfile1)
    save(auc, file=outfile2)
    save(bin, file=outfile3)
    save(f, file=outfile4)
}

BIN <- function(out, x){
	L <- scTensor:::.HCLUST(out$ligand[x[1],])
	R <- scTensor:::.HCLUST(out$receptor[x[2],])
	LR <- scTensor:::.HCLUST(out$lrpair[x[3],])
    # L <- .DBSCAN(out$ligand[x[1],])
    # R <- .DBSCAN(out$receptor[x[2],])
    # LR <- .DBSCAN(out$lrpair[x[3],])
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

ROC_AUC_BIN_F_previous_sctensor <- function(trueCaH, out, outfile1, outfile2, outfile3, outfile4){
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
    rocaucf <- ROCAUCF(trueCaH, score, bin)
    roc <- lapply(rocaucf, function(raf){
    	raf$roc
    })
    auc <- lapply(rocaucf, function(raf){
    	raf$auc
    })
    f <- lapply(rocaucf, function(raf){
    	raf$f
    })
    # Save
    save(roc, file=outfile1)
    save(auc, file=outfile2)
    save(bin, file=outfile3)
    save(f, file=outfile4)
}

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
    L <- scTensor:::.HCLUST(out$ligand[x[1],])
    R <- scTensor:::.HCLUST(out$receptor[x[2],])
    L[which(L == "selected")] <- 1
    R[which(R == "selected")] <- 1
    L[which(L == "not selected")] <- 0
    R[which(R == "not selected")] <- 0
    L <- as.numeric(L)
    R <- as.numeric(R)
    list(ligand=L, receptor=R)
}

ROC_AUC_BIN_F_sctensor <- function(trueCaH, out, outfile1, outfile2, outfile3, outfile4){
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
        tmp <- scTensor:::.HCLUST(as.vector(bin.core[,,i]))
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
    rocaucf <- ROCAUCF(trueCaH, score, bin)
    roc <- lapply(rocaucf, function(raf){
    	raf$roc
    })
    auc <- lapply(rocaucf, function(raf){
    	raf$auc
    })
    f <- lapply(rocaucf, function(raf){
    	raf$f
    })
    # Save
    save(roc, file=outfile1)
    save(auc, file=outfile2)
    save(bin, file=outfile3)
    save(f, file=outfile4)
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
