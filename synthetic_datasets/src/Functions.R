library("scTensor")
library("scales")
library("RColorBrewer")
library("ROCR")
library("nnTensor")
library("rTensor")
library("ggplot2")
library("viridis")
library("dbscan")

# Time
Plot_Memory <- function(method, E, outfile){
    df <- aggregateMemory(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=GB))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="GB", limits=c(0, 15))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=10, height=10)
}

aggregateMemory <- function(E, Method){
    out <- c()
    Celltypes=c("3","5","10","20","30")
    CCIPatterns=c("1","3","5")
    CCI_Complexity=c("OnetoOne", "ManytoMany")
    nCL <- length(Celltypes)
    nCP <- length(CCIPatterns)
    nCC <- length(CCI_Complexity)
    for(i in 1:nCL){
        for(j in 1:nCP){
            for(k in 1:nCC){
                if(Method %in% c("sctensor", "previous_sctensor")){
                    inputfile <- paste0("benchmarks/",Method, "_", E, "_1_",
                        Celltypes[i], "Celltypes_",
                        CCIPatterns[j], "CCIPatterns_",
                        CCI_Complexity[k], ".txt")
                }else{
                    inputfile <- paste0("benchmarks/",
                        E, "_",                         
                        Method, "_",
                        Celltypes[i], "Celltypes_",
                        CCIPatterns[j], "CCIPatterns_",
                        CCI_Complexity[k], ".txt")      
                }
                memory <- read.delim(inputfile)
                out <- c(out, memory$max_rss / 10^3)
            }
        }
    }
    tmp1 <- as.vector(
            sapply(Celltypes, function(x){
                rep(x, length=nCP*nCC)}))
    tmp2 <- as.vector(
            unlist(sapply(Celltypes, function(y){
                unlist(sapply(paste0(CCIPatterns, "CCI_Patterns"), function(x){
                    rep(x, length=nCC)
                }))})))
    tmp2 <- gsub("CCI_Patterns", "", tmp2)
    tmp3 <- rep(CCI_Complexity, length=nCL*nCP*nCC)
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        GB=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df
}

# Time
Plot_Time <- function(method, E, outfile){
    df <- aggregateTime(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=Time))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="Hour", limits=c(0, 55))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=10, height=10)
}

aggregateTime <- function(E, Method){
    out <- c()
    Celltypes=c("3","5","10","20","30")
    CCIPatterns=c("1","3","5")
    CCI_Complexity=c("OnetoOne", "ManytoMany")
    nCL <- length(Celltypes)
    nCP <- length(CCIPatterns)
    nCC <- length(CCI_Complexity)
    for(i in 1:nCL){
        for(j in 1:nCP){
            for(k in 1:nCC){
                if(Method %in% c("sctensor", "previous_sctensor")){
                    inputfile <- paste0("benchmarks/",
                        Method, "_", E, "_1_",
                        Celltypes[i], "Celltypes_",
                        CCIPatterns[j], "CCIPatterns_",
                        CCI_Complexity[k], ".txt")
                }else{
                    inputfile <- paste0("benchmarks/",
                        E, "_", 
                        Method, "_",
                        Celltypes[i], "Celltypes_",
                        CCIPatterns[j], "CCIPatterns_",
                        CCI_Complexity[k], ".txt")      
                }
                time <- read.delim(inputfile)
                out <- c(out, time$s / 60 / 60)
            }
        }
    }
    tmp1 <- as.vector(
            sapply(Celltypes, function(x){
                rep(x, length=nCP*nCC)}))
    tmp2 <- as.vector(
            unlist(sapply(Celltypes, function(y){
                unlist(sapply(paste0(CCIPatterns, "CCI_Patterns"), function(x){
                    rep(x, length=nCC)
                }))})))
    tmp2 <- gsub("CCI_Patterns", "", tmp2)
    tmp3 <- rep(CCI_Complexity, length=nCL*nCP*nCC)
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        Time=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df
}

# F
Plot_F <- function(method, E, outfile){
    df <- aggregateF(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=F))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="F-measure", limits=c(0, 0.7))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=10, height=10)
}

aggregateF <- function(E, Method){
    out <- c()
    Celltypes=c("3","5","10","20","30")
    CCIPatterns=c("1","3","5")
    CCI_Complexity=c("OnetoOne", "ManytoMany")
    nCL <- length(Celltypes)
    nCP <- length(CCIPatterns)
    nCC <- length(CCI_Complexity)
    for(i in 1:nCL){
        for(j in 1:nCP){
            for(k in 1:nCC){
                inputfile <- paste0("output/",
                    E, "/",
                    Method, "/F/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                f <- unlist(f)
                target <- which(!is.nan(f))
                out = c(out, mean(f[target]))
            }
        }
    }
    tmp1 <- as.vector(
            sapply(Celltypes, function(x){
                rep(x, length=nCP*nCC)}))
    tmp2 <- as.vector(
            unlist(sapply(Celltypes, function(y){
                unlist(sapply(paste0(CCIPatterns, "CCI_Patterns"), function(x){
                    rep(x, length=nCC)
                }))})))
    tmp2 <- gsub("CCI_Patterns", "", tmp2)
    tmp3 <- rep(CCI_Complexity, length=nCL*nCP*nCC)
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        F=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df
}

# AUC
Plot_AUC <- function(method, E, outfile){
    df <- aggregateAUC(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=AUC))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="AUC", limits=c(0, 1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=10, height=10)
}

aggregateAUC <- function(E, Method){
    out <- c()
    Celltypes=c("3","5","10","20","30")
    CCIPatterns=c("1","3","5")
    CCI_Complexity=c("OnetoOne", "ManytoMany")
    nCL <- length(Celltypes)
    nCP <- length(CCIPatterns)
    nCC <- length(CCI_Complexity)
    for(i in 1:nCL){
        for(j in 1:nCP){
            for(k in 1:nCC){
                inputfile <- paste0("output/",
                    E, "/",
                    Method, "/AUC/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                auc <- unlist(auc)
                target <- which(!is.nan(auc))
                out = c(out, mean(auc[target]))
            }
        }
    }
    tmp1 <- as.vector(
            sapply(Celltypes, function(x){
                rep(x, length=nCP*nCC)}))
    tmp2 <- as.vector(
            unlist(sapply(Celltypes, function(y){
                unlist(sapply(paste0(CCIPatterns, "CCI_Patterns"), function(x){
                    rep(x, length=nCC)
                }))})))
    tmp2 <- gsub("CCI_Patterns", "", tmp2)
    tmp3 <- rep(CCI_Complexity, length=nCL*nCP*nCC)
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        AUC=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df
}

Plot_ROC_AUC_F <- function(roc, auc, f, outfile){
    colvec <- c(brewer.pal(9, "Set1"),
        brewer.pal(8, "Set2"), 
        brewer.pal(8, "Dark2"))
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
    if(11 <= l && l < 25){
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
	cif <- cciInfo[["E2"]][[spl]]
	target <- length(grep("CCI", names(cif)))
	ncelltypes <- length(cif$CCI1$LPattern)
	trueCaH <- list()
	length(trueCaH) <- target
	counter <- 1
	# Each CCI1, 2, 3, ...
	for(j in seq(target)){
	    L <- eval(parse(text=paste0("cif$CCI", j, "$LPattern")))
	    R <- eval(parse(text=paste0("cif$CCI", j, "$RPattern")))
	    nGene <- eval(parse(text=paste0("cif$CCI", j, "$nGene")))
	    # trueCaH
	    tmpTrueCaH <- rep(as.vector(base::outer(L, R, "*")), nGene)
	    # Setting true CaH as 1
	    trueCaH[[j]] <- rep(0, length=ncelltypes^2*cif$nPair)
	    trueCaH[[j]][counter:(counter+length(tmpTrueCaH)-1)] <- tmpTrueCaH

	    # Setting color scheme
	    tmpTrueColor <- ifelse(tmpTrueCaH == 0, alpha(rgb(0,0,0), 0.1), brewer.pal(9, "Set1")[j])
	    counter <- counter + length(tmpTrueCaH)
	}
	list(trueCaH=trueCaH, ncelltypes=ncelltypes, cif=cif, spl=spl)
}

cciInfo <- list(
    "E2"=list(
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
    ),
    "E5"=list(
        ###### 1. 3 Celltypes, 1 CCI Patterns, One-by-one
        "3Celltypes_1CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,0),
             nGene=50,
             fc="E5")
          ),
        ###### 2. 5 Celltypes, 1 CCI Patterns, One-by-one
        "5Celltypes_1CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0),
             RPattern=c(0,1,0,0,0),
             nGene=50,
             fc="E5")
          ),
        ###### 3. 10 Celltypes, 1 CCI Patterns, One-by-one
        "10Celltypes_1CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5")
          ),
        ###### 4. 20 Celltypes, 1 CCI Patterns, One-by-one
        "20Celltypes_1CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5")
          ),
        ###### 5. 30 Celltypes, 1 CCI Patterns, One-by-one
        "30Celltypes_1CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5")
          ),

        ###### 6. 3 Celltypes, 3 CCI Patterns, One-by-one
        # Fig.のCase I
        "3Celltypes_3CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,1,0),
             RPattern=c(0,0,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(0,0,1),
             RPattern=c(1,0,0),
             nGene=50,
             fc="E5")
          ),
        ###### 7. 5 Celltypes, 3 CCI Patterns, One-by-one
        "5Celltypes_3CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0),
             RPattern=c(0,1,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,1,0,0,0),
             RPattern=c(0,0,1,0,0),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(0,0,1,0,0),
             RPattern=c(0,0,0,1,0),
             nGene=50,
             fc="E5")
          ),
        ###### 8. 10 Celltypes, 3 CCI Patterns, One-by-one
        "10Celltypes_3CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0),
             nGene=50,
             fc="E5")
          ),
        ###### 9. 20 Celltypes, 3 CCI Patterns, One-by-one
        "20Celltypes_3CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5")
          ),
        ###### 10. 30 Celltypes, 3 CCI Patterns, One-by-one
        "30Celltypes_3CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5")
          ),

        ###### 11. 3 Celltypes, 5 CCI Patterns, One-by-one
        "3Celltypes_5CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,1,0),
             RPattern=c(0,0,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(0,0,1),
             RPattern=c(1,0,0),
             nGene=50,
             fc="E5"),
          CCI4=list(
             LPattern=c(0,1,0),
             RPattern=c(1,0,0),
             nGene=50,
             fc="E5"),
          CCI5=list(
             LPattern=c(0,0,1),
             RPattern=c(0,1,0),
             nGene=50,
             fc="E5")
          ),
        ###### 12. 5 Celltypes, 5 CCI Patterns, One-by-one
        "5Celltypes_5CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0),
             RPattern=c(0,1,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,1,0,0,0),
             RPattern=c(0,0,1,0,0),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(0,0,1,0,0),
             RPattern=c(0,0,0,1,0),
             nGene=50,
             fc="E5"),
          CCI4=list(
             LPattern=c(0,0,0,1,0),
             RPattern=c(0,0,0,0,1),
             nGene=50,
             fc="E5"),
          CCI5=list(
             LPattern=c(0,0,0,0,1),
             RPattern=c(1,0,0,0,0),
             nGene=50,
             fc="E5")
          ),
        ###### 13. 10 Celltypes, 5 CCI Patterns, One-by-one
        "10Celltypes_5CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI4=list(
             LPattern=c(0,0,0,1,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,1,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI5=list(
             LPattern=c(0,0,0,0,1,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,0,0,0,0),
             nGene=50,
             fc="E5")
          ),
        ###### 14. 20 Celltypes, 5 CCI Patterns, One-by-one
        "20Celltypes_5CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI4=list(
             LPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI5=list(
             LPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5")
          ),
        ###### 15. 30 Celltypes, 5 CCI Patterns, One-by-one
        "30Celltypes_5CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI4=list(
             LPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI5=list(
             LPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5")
          ),

        ###### 16. 3 Celltypes, 1 CCI Patterns, Many-by-many
        "3Celltypes_1CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,1),
             nGene=50,
             fc="E5")
          ),
        ###### 17. 5 Celltypes, 1 CCI Patterns, Many-by-many
        "5Celltypes_1CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,0,0,0),
             RPattern=c(0,0,1,1,1),
             nGene=50,
             fc="E5")
          ),
        ###### 18. 10 Celltypes, 1 CCI Patterns, Many-by-many
        "10Celltypes_1CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1),
             nGene=50,
             fc="E5")
          ),
        ###### 19. 20 Celltypes, 1 CCI Patterns, Many-by-many
        "20Celltypes_1CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5")
          ),
        ###### 20. 30 Celltypes, 1 CCI Patterns, Many-by-many
        "30Celltypes_1CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5")
          ),

        ###### 21. 3 Celltypes, 3 CCI Patterns, Many-by-many
        # Fig.のCase II
        "3Celltypes_3CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,1),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,0,1),
             RPattern=c(1,1,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(1,1,1),
             RPattern=c(1,0,1),
             nGene=50,
             fc="E5")
          ),
        ###### 22. 5 Celltypes, 3 CCI Patterns, Many-by-many
        "5Celltypes_3CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,0,0,0),
             RPattern=c(0,0,1,1,1),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,0,0,0,1),
             RPattern=c(1,1,1,1,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(1,1,1,1,1),
             RPattern=c(1,1,0,0,1),
             nGene=50,
             fc="E5")
          ),
        ###### 23. 10 Celltypes, 3 CCI Patterns, Many-by-many
        "10Celltypes_3CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E5")
          ),
        ###### 24. 20 Celltypes, 3 CCI Patterns, Many-by-many
        "20Celltypes_3CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E5")
          ),
        ###### 25. 30 Celltypes, 3 CCI Patterns, Many-by-many
        "30Celltypes_3CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E5")
          ),

        ###### 26. 3 Celltypes, 5 CCI Patterns, Many-by-many
        "3Celltypes_5CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,1),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,0,1),
             RPattern=c(1,1,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(1,1,1),
             RPattern=c(1,0,1),
             nGene=50,
             fc="E5"),
          CCI4=list(
             LPattern=c(1,0,1),
             RPattern=c(0,1,0),
             nGene=50,
             fc="E5"),
          CCI5=list(
             LPattern=c(0,1,0),
             RPattern=c(1,0,1),
             nGene=50,
             fc="E5")
          ),
        ###### 27. 5 Celltypes, 5 CCI Patterns, Many-by-many
        "5Celltypes_5CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,0,0,0),
             RPattern=c(0,0,1,1,1),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,0,0,0,1),
             RPattern=c(1,1,1,1,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(1,1,1,1,1),
             RPattern=c(1,1,0,0,1),
             nGene=50,
             fc="E5"),
          CCI4=list(
             LPattern=c(1,0,1,0,1),
             RPattern=c(0,1,0,1,0),
             nGene=50,
             fc="E5"),
          CCI5=list(
             LPattern=c(0,1,0,1,0),
             RPattern=c(1,0,1,0,1),
             nGene=50,
             fc="E5")
          ),
        ###### 28. 10 Celltypes, 5 CCI Patterns, Many-by-many
        "10Celltypes_5CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E5"),
          CCI4=list(
             LPattern=c(1,0,1,0,1,0,1,0,1,0),
             RPattern=c(0,1,0,1,0,1,0,1,0,1),
             nGene=50,
             fc="E5"),
          CCI5=list(
             LPattern=c(0,1,0,1,0,1,0,1,0,1),
             RPattern=c(1,0,1,0,1,0,1,0,1,0),
             nGene=50,
             fc="E5")
          ),
        ###### 29. 20 Celltypes, 5 CCI Patterns, Many-by-many
        "20Celltypes_5CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E5"),
          CCI4=list(
             LPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
             RPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
             nGene=50,
             fc="E5"),
          CCI5=list(
             LPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
             RPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
             nGene=50,
             fc="E5")
          ),
        ###### 30. 30 Celltypes, 5 CCI Patterns, Many-by-many
        "30Celltypes_5CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E5"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E5"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E5"),
          CCI4=list(
             LPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
             RPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
             nGene=50,
             fc="E5"),
          CCI5=list(
             LPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
             RPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
             nGene=50,
             fc="E5")
          )
    ),
    "E10"=list(
        ###### 1. 3 Celltypes, 1 CCI Patterns, One-by-one
        "3Celltypes_1CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,0),
             nGene=50,
             fc="E10")
          ),
        ###### 2. 5 Celltypes, 1 CCI Patterns, One-by-one
        "5Celltypes_1CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0),
             RPattern=c(0,1,0,0,0),
             nGene=50,
             fc="E10")
          ),
        ###### 3. 10 Celltypes, 1 CCI Patterns, One-by-one
        "10Celltypes_1CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10")
          ),
        ###### 4. 20 Celltypes, 1 CCI Patterns, One-by-one
        "20Celltypes_1CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10")
          ),
        ###### 5. 30 Celltypes, 1 CCI Patterns, One-by-one
        "30Celltypes_1CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10")
          ),

        ###### 6. 3 Celltypes, 3 CCI Patterns, One-by-one
        # Fig.のCase I
        "3Celltypes_3CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,1,0),
             RPattern=c(0,0,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(0,0,1),
             RPattern=c(1,0,0),
             nGene=50,
             fc="E10")
          ),
        ###### 7. 5 Celltypes, 3 CCI Patterns, One-by-one
        "5Celltypes_3CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0),
             RPattern=c(0,1,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,1,0,0,0),
             RPattern=c(0,0,1,0,0),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(0,0,1,0,0),
             RPattern=c(0,0,0,1,0),
             nGene=50,
             fc="E10")
          ),
        ###### 8. 10 Celltypes, 3 CCI Patterns, One-by-one
        "10Celltypes_3CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0),
             nGene=50,
             fc="E10")
          ),
        ###### 9. 20 Celltypes, 3 CCI Patterns, One-by-one
        "20Celltypes_3CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10")
          ),
        ###### 10. 30 Celltypes, 3 CCI Patterns, One-by-one
        "30Celltypes_3CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10")
          ),

        ###### 11. 3 Celltypes, 5 CCI Patterns, One-by-one
        "3Celltypes_5CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,1,0),
             RPattern=c(0,0,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(0,0,1),
             RPattern=c(1,0,0),
             nGene=50,
             fc="E10"),
          CCI4=list(
             LPattern=c(0,1,0),
             RPattern=c(1,0,0),
             nGene=50,
             fc="E10"),
          CCI5=list(
             LPattern=c(0,0,1),
             RPattern=c(0,1,0),
             nGene=50,
             fc="E10")
          ),
        ###### 12. 5 Celltypes, 5 CCI Patterns, One-by-one
        "5Celltypes_5CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0),
             RPattern=c(0,1,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,1,0,0,0),
             RPattern=c(0,0,1,0,0),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(0,0,1,0,0),
             RPattern=c(0,0,0,1,0),
             nGene=50,
             fc="E10"),
          CCI4=list(
             LPattern=c(0,0,0,1,0),
             RPattern=c(0,0,0,0,1),
             nGene=50,
             fc="E10"),
          CCI5=list(
             LPattern=c(0,0,0,0,1),
             RPattern=c(1,0,0,0,0),
             nGene=50,
             fc="E10")
          ),
        ###### 13. 10 Celltypes, 5 CCI Patterns, One-by-one
        "10Celltypes_5CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI4=list(
             LPattern=c(0,0,0,1,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,1,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI5=list(
             LPattern=c(0,0,0,0,1,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,0,0,0,0),
             nGene=50,
             fc="E10")
          ),
        ###### 14. 20 Celltypes, 5 CCI Patterns, One-by-one
        "20Celltypes_5CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI4=list(
             LPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI5=list(
             LPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10")
          ),
        ###### 15. 30 Celltypes, 5 CCI Patterns, One-by-one
        "30Celltypes_5CCIPatterns_OnetoOne"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI4=list(
             LPattern=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI5=list(
             LPattern=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10")
          ),

        ###### 16. 3 Celltypes, 1 CCI Patterns, Many-by-many
        "3Celltypes_1CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,1),
             nGene=50,
             fc="E10")
          ),
        ###### 17. 5 Celltypes, 1 CCI Patterns, Many-by-many
        "5Celltypes_1CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,0,0,0),
             RPattern=c(0,0,1,1,1),
             nGene=50,
             fc="E10")
          ),
        ###### 18. 10 Celltypes, 1 CCI Patterns, Many-by-many
        "10Celltypes_1CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1),
             nGene=50,
             fc="E10")
          ),
        ###### 19. 20 Celltypes, 1 CCI Patterns, Many-by-many
        "20Celltypes_1CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10")
          ),
        ###### 20. 30 Celltypes, 1 CCI Patterns, Many-by-many
        "30Celltypes_1CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10")
          ),

        ###### 21. 3 Celltypes, 3 CCI Patterns, Many-by-many
        # Fig.のCase II
        "3Celltypes_3CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,1),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,0,1),
             RPattern=c(1,1,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(1,1,1),
             RPattern=c(1,0,1),
             nGene=50,
             fc="E10")
          ),
        ###### 22. 5 Celltypes, 3 CCI Patterns, Many-by-many
        "5Celltypes_3CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,0,0,0),
             RPattern=c(0,0,1,1,1),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,0,0,0,1),
             RPattern=c(1,1,1,1,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(1,1,1,1,1),
             RPattern=c(1,1,0,0,1),
             nGene=50,
             fc="E10")
          ),
        ###### 23. 10 Celltypes, 3 CCI Patterns, Many-by-many
        "10Celltypes_3CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E10")
          ),
        ###### 24. 20 Celltypes, 3 CCI Patterns, Many-by-many
        "20Celltypes_3CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E10")
          ),
        ###### 25. 30 Celltypes, 3 CCI Patterns, Many-by-many
        "30Celltypes_3CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E10")
          ),

        ###### 26. 3 Celltypes, 5 CCI Patterns, Many-by-many
        "3Celltypes_5CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,0,0),
             RPattern=c(0,1,1),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,0,1),
             RPattern=c(1,1,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(1,1,1),
             RPattern=c(1,0,1),
             nGene=50,
             fc="E10"),
          CCI4=list(
             LPattern=c(1,0,1),
             RPattern=c(0,1,0),
             nGene=50,
             fc="E10"),
          CCI5=list(
             LPattern=c(0,1,0),
             RPattern=c(1,0,1),
             nGene=50,
             fc="E10")
          ),
        ###### 27. 5 Celltypes, 5 CCI Patterns, Many-by-many
        "5Celltypes_5CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,0,0,0),
             RPattern=c(0,0,1,1,1),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,0,0,0,1),
             RPattern=c(1,1,1,1,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(1,1,1,1,1),
             RPattern=c(1,1,0,0,1),
             nGene=50,
             fc="E10"),
          CCI4=list(
             LPattern=c(1,0,1,0,1),
             RPattern=c(0,1,0,1,0),
             nGene=50,
             fc="E10"),
          CCI5=list(
             LPattern=c(0,1,0,1,0),
             RPattern=c(1,0,1,0,1),
             nGene=50,
             fc="E10")
          ),
        ###### 28. 10 Celltypes, 5 CCI Patterns, Many-by-many
        "10Celltypes_5CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E10"),
          CCI4=list(
             LPattern=c(1,0,1,0,1,0,1,0,1,0),
             RPattern=c(0,1,0,1,0,1,0,1,0,1),
             nGene=50,
             fc="E10"),
          CCI5=list(
             LPattern=c(0,1,0,1,0,1,0,1,0,1),
             RPattern=c(1,0,1,0,1,0,1,0,1,0),
             nGene=50,
             fc="E10")
          ),
        ###### 29. 20 Celltypes, 5 CCI Patterns, Many-by-many
        "20Celltypes_5CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E10"),
          CCI4=list(
             LPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
             RPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
             nGene=50,
             fc="E10"),
          CCI5=list(
             LPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
             RPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
             nGene=50,
             fc="E10")
          ),
        ###### 30. 30 Celltypes, 5 CCI Patterns, Many-by-many
        "30Celltypes_5CCIPatterns_ManytoMany"=
        list(nPair=1000,
          CCI1=list(
             LPattern=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             RPattern=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             nGene=50,
             fc="E10"),
          CCI2=list(
             LPattern=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
             RPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             nGene=50,
             fc="E10"),
          CCI3=list(
             LPattern=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             RPattern=c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
             nGene=50,
             fc="E10"),
          CCI4=list(
             LPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
             RPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
             nGene=50,
             fc="E10"),
          CCI5=list(
             LPattern=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),
             RPattern=c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),
             nGene=50,
             fc="E10")
          )
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
