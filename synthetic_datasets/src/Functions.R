library("igraph")
library("scTensor")
library("nnTensor")
library("ccTensor")
library("scales")
library("RColorBrewer")
library("ROCR")
library("rTensor")
library("ggplot2")
library("viridis")
library("dbscan")
library("fields")
library("mclust")

# Check 2.4.0
if(packageVersion("scTensor") != "2.4.0"){
    stop("scTensor is not 2.4.0!")
}

# Plot Ground Truth
mylayout <- function(l){
    rbind(
        cbind(2, rev(seq(l))),
        cbind(1, rev(seq(l)))
    )
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
    myedgecolor <- c(rgb(1,0,0), rgb(0,1,0),
        rgb(0,0,1), rgb(1,0,1), rgb(0.5,0,1))
    plot(g, layout = mylayout(l),
        vertex.color=mynodecolor[V(g)$type + 1],
        edge.color=myedgecolor[x-1])
}

# L-R pairs
# Method <- "halpern"
# E <- "E10"
Plot_LR <- function(Method, E, outfile){
    Celltypes=c("3","5","10","20","30")
    CCIPatterns=c("1","3","5")
    CCI_Complexity=c("OnetoOne", "ManytoMany")
    nCL <- length(Celltypes)
    nCP <- length(CCIPatterns)
    nCC <- length(CCI_Complexity)
    for(i in 1:nCL){
        for(j in 1:nCP){
            for(k in 1:nCC){
                # Data loading
                # out, 1-out@pval, out@tnsr
                inputfile1 <- paste0("output/",
                    E, "/",
                    Method, "/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                # Binarized
                inputfile2 <- paste0("output/",
                    E, "/",
                    Method, "/BIN/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                # cif
                inputfile3 <- paste0("data/groundtruth/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile1)
                load(inputfile2)
                load(inputfile3)
                plotSlice(out, bin, trueCaH,
                    ncelltypes, E, Method, cif,
                    Celltypes[i], CCIPatterns[j], CCI_Complexity[k])
            }
        }
    }
    file.create(outfile)
}

# celltypes=Celltypes[i]
# ccipatterns=CCIPatterns[j]
# cci_complexity=CCI_Complexity[k]
plotSlice <- function(out, bin, trueCaH, ncelltypes, E, Method, cif, celltypes, ccipatterns, cci_complexity){
    tensor.methods <- c(
        "labelpermutation_tensor", "labelpermutation2_tensor",
        "halpern_tensor", "cabelloaguilar_tensor")
    pval.methods <- c(
        "labelpermutation", "labelpermutation2",
        "halpern", "cabelloaguilar")
    if(Method %in% tensor.methods){
        g <- out$tnsr
        bintnsr <- Vec2Tensor(bin, ncelltypes)
    }
    if(Method %in% pval.methods){
        g <- 1 - out$pval
        bintnsr <- Vec2Tensor(bin, ncelltypes)
    }
    if(Method == "previous_sctensor"){
        g <- returnTensor(out, cif, Method)@data
        # Scoring
        score <- apply(out$index, 1, function(x){
            nnTensor::recTensor(x[4],
                list(as.matrix(out$ligand[x[1], ]),
                    as.matrix(out$receptor[x[2], ]),
                    as.matrix(out$lrpair[x[3], ])),
                        reverse=TRUE)
        })
    }
    if(Method == "sctensor"){
        g <- returnTensor(out, cif, Method)@data
        # Scoring
        score <- apply(out$index, 1, function(x){
            nnTensor::recTensor(x[4],
                list(as.matrix(out$ligand[x[1], ]),
                    as.matrix(out$receptor[x[2], ]),
                    as.matrix(out$lrpair[x[1], x[2], ]@data)), reverse=TRUE)
        })
    }
    # Plot
    counter <- 1 # 1,51,101,...
    for(l in 2:length(cif)){
        # Setting
        outdir <- paste0("plot/",
            E, "/L-R/",
            Method, "/",
            celltypes, "Celltypes_",
            ccipatterns, "CCIPatterns_",
            cci_complexity, "/",
            "CCI", l-1)
        dir.create(outdir, recursive=TRUE)
        # 1. True CCI
        truematrix <- outer(cif[[l]]$LPattern, cif[[l]]$RPattern)
        truematrix <- t(truematrix[nrow(truematrix):1,])
        match_tnsr_method <- grep("_tensor", Method)
        if(length(match_tnsr_method) == 0){
            # 2. Modesum of score tensor
            start <- counter
            end <- counter + cif[[l]]$nGene - 1
            outmatrix1 <- modeSum(
                as.tensor(g[,,start:end]), m=3, drop=TRUE)@data
            outmatrix1 <- t(outmatrix1[nrow(outmatrix1):1,])
            # 3. Modesum of bin tensor
            if(Method %in% c("previous_sctensor", "sctensor")){
                bintnsr <- Vec2Tensor_sctensor(score, bin, trueCaH[[l-1]], ncelltypes)
            }
            outmatrix2 <- modeSum(
                as.tensor(bintnsr[,,start:end]), m=3, drop=TRUE)@data
            outmatrix2 <- t(outmatrix2[nrow(outmatrix2):1,])
            # Plot
            outfile1 <- paste0(outdir, "/CCI", l-1, ".png")
            png(file=outfile1, width=2100, height=700)
            layout(t(1:3))
            image.plot(outmatrix1, xlab="R", ylab="L",
                main="Score")
            image.plot(outmatrix2, xlab="R", ylab="L",
                main="Binarization")
            image.plot(truematrix, xlab="R", ylab="L",
                col=c(rgb(0,0,0,0), rgb(1,0,1)),
                main="True CCI")
            dev.off()
        }
        counter2 <- counter # 1,2,3,...51,52,53,...
        counter <- counter + cif[[l]]$nGene
        # Each L-R pair
        for(m in seq(cif[[l]]$nGene)){
            outfile2 <- paste0(outdir, "/", counter2, ".png")
            outmatrix3 <- t(g[nrow(g):1,, counter2])
            # Plot
            png(file=outfile2, width=1400, height=700)
            layout(t(1:2))
            image.plot(outmatrix3, xlab="R", ylab="L",
                main="Score")
            image.plot(truematrix, xlab="R", ylab="L",
                col=c(rgb(0,0,0,0), rgb(1,0,1)),
                main="True CCI")
           dev.off()
            counter2 <- counter2 + 1
        }
    }
}

returnTensor <- function(out, cif, Method){
    g <- array(0,
        dim=c(length(cif$CCI1$LPattern),
            length(cif$CCI1$RPattern), cif$nPair))
    if(Method == "previous_sctensor"){
        for(i in seq_len(nrow(out$index))){
            mode1 <- out$index[i, "Mode1"]
            mode2 <- out$index[i, "Mode2"]
            mode3 <- out$index[i, "Mode3"]
            core <- out$index[i, "Value"]
            g <- g + nnTensor::recTensor(core,
                list(as.matrix(out$ligand[mode1, ]),
                    as.matrix(out$receptor[mode2, ]),
                    as.matrix(out$lrpair[mode3, ])),
                    reverse=TRUE)
        }
    }
    if(Method == "sctensor"){
        for(i in seq_len(nrow(out$index))){
            mode1 <- out$index[i, "Mode1"]
            mode2 <- out$index[i, "Mode2"]
            core <- out$index[i, "Value"]
            g <- g + nnTensor::recTensor(core,
                list(as.matrix(out$ligand[mode1, ]),
                    as.matrix(out$receptor[mode2, ]),
                    as.matrix(out$lrpair[mode1, mode2, ]@data)),
                    reverse=TRUE)
        }
    }
    g
}

# Memory
Plot_Memory <- function(method, E, outfile){
    df <- aggregateMemory(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=GB))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="GB", limits=c(0, 17.5))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
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
                out <- c(out, memory$max_rss)
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

aggregateMemory_merge <- function(E){
    Methods <- c(
        "labelpermutation_tensor",
        "labelpermutation",
        "labelpermutation2_tensor",
        "labelpermutation2",
        "halpern_tensor",
        "halpern",
        "cabelloaguilar_tensor",
        "cabelloaguilar",
        "previous_sctensor",
        "sctensor")
    df <- data.frame()
    for(i in seq_along(Methods)){
        tmp.df <- aggregateMemory(E, Methods[i])
        tmp.df <- cbind(tmp.df, Methods[i])
        df <- rbind(df, tmp.df)
    }
    colnames(df)[5] <- "Methods"
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Methods <- factor(df$Methods, levels=Methods)
    df
}

Plot_Memory_Merge <- function(E, outfile){
    df <- aggregateMemory_merge(E=E)
    # Plot
    gg <- ggplot(df, aes(x=Methods, y=GB, fill=Methods))
    gg <- gg + geom_boxplot()
    gg <- gg + facet_wrap(~CCI_Complexity)
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.title.x=element_blank())
    gg <- gg + theme(axis.text.x=element_blank())
    gg <- gg + theme(axis.title.y=element_text(size = 30))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(1, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=20, height=10)
}

# Time
Plot_Time <- function(method, E, outfile){
    df <- aggregateTime(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=Time))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="Hour", limits=c(0, 45))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
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

aggregateTime_merge <- function(E){
    Methods <- c(
        "labelpermutation_tensor",
        "labelpermutation",
        "labelpermutation2_tensor",
        "labelpermutation2",
        "halpern_tensor",
        "halpern",
        "cabelloaguilar_tensor",
        "cabelloaguilar",
        "previous_sctensor",
        "sctensor")
    df <- data.frame()
    for(i in seq_along(Methods)){
        tmp.df <- aggregateTime(E, Methods[i])
        tmp.df <- cbind(tmp.df, Methods[i])
        df <- rbind(df, tmp.df)
    }
    colnames(df)[5] <- "Methods"
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Methods <- factor(df$Methods, levels=Methods)
    df
}

Plot_Time_Merge <- function(E, outfile){
    df <- aggregateTime_merge(E=E)
    # Plot
    gg <- ggplot(df, aes(x=Methods, y=Time, fill=Methods))
    gg <- gg + geom_boxplot()
    gg <- gg + facet_wrap(~CCI_Complexity)
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.title.x=element_blank())
    gg <- gg + theme(axis.text.x=element_blank())
    gg <- gg + theme(axis.title.y=element_text(size = 30))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(1, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=20, height=10)
}


# F
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
    df$F[which(is.nan(df$F))] <- 0
    df$F[which(is.na(df$F))] <- 0
    df
}

aggregateF_eachCCI <- function(E, Method){
    out <- c()
    tmp1 <- c()
    tmp2 <- c()
    tmp3 <- c()
    tmp4 <- c()
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
                out = c(out, f)
                l <- length(f)
                tmp1 <- c(tmp1, rep(Celltypes[i], l))
                tmp2 <- c(tmp2, rep(CCIPatterns[j], l))
                tmp3 <- c(tmp3, rep(CCI_Complexity[k], l))
                tmp4 <- c(tmp4, seq(l))
            }
        }
    }
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        Replicates=tmp4,
        F=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$Replicates <- as.character(df$Replicates)
    df$F[which(is.nan(df$F))] <- 0
    df$F[which(is.na(df$F))] <- 0
    df
}

Plot_F <- function(method, E, outfile){
    df <- aggregateF(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=F))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="F-measure", limits=c(0, 1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

Plot_F_eachCCI <- function(method, E, outfile){
    df <- aggregateF_eachCCI(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=Replicates, fill=F))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~CCI_Complexity, scales="free_y", space="free_y")
    gg <- gg + scale_fill_viridis(name="F-measure", limits=c(0, 1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

aggregateF_merge <- function(E){
    BinMethods <- c(
        "labelpermutation",
        "labelpermutation2",
        "halpern",
        "cabelloaguilar",
        "previous_sctensor",
        "sctensor")
    df <- data.frame()
    for(i in seq_along(BinMethods)){
        tmp.df <- aggregateF_eachCCI(E, BinMethods[i])
        tmp.df <- cbind(tmp.df, BinMethods[i])
        df <- rbind(df, tmp.df)
    }
    colnames(df)[6] <- "BinMethods"
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$BinMethods <- factor(df$BinMethods, levels=BinMethods)
    df
}

Plot_F_Merge <- function(E, outfile){
    df <- aggregateF_merge(E=E)
    # Plot
    gg <- ggplot(df, aes(x=BinMethods, y=F, fill=BinMethods))
    gg <- gg + geom_boxplot()
    gg <- gg + facet_wrap(~CCI_Complexity)
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.title.x=element_blank())
    gg <- gg + theme(axis.text.x=element_blank())
    gg <- gg + theme(axis.title.y=element_text(size = 30))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(1, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=20, height=10)
}


# MCC
aggregateMCC <- function(E, Method){
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
                    Method, "/MCC/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                mcc <- unlist(mcc)
                target <- which(!is.nan(mcc))
                out = c(out, mean(mcc[target]))
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
        MCC=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$MCC[which(is.nan(df$MCC))] <- NA
    df$MCC[which(df$MCC < 0)] <- 0
    df
}

aggregateMCC_eachCCI <- function(E, Method){
    out <- c()
    tmp1 <- c()
    tmp2 <- c()
    tmp3 <- c()
    tmp4 <- c()
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
                    Method, "/MCC/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                mcc <- unlist(mcc)
                target <- which(!is.nan(mcc))
                l <- length(target)
                out = c(out, mcc[target])
                tmp1 <- c(tmp1, rep(Celltypes[i], l))
                tmp2 <- c(tmp2, rep(CCIPatterns[j], l))
                tmp3 <- c(tmp3, rep(CCI_Complexity[k], l))
                if(l != 0){
                    tmp4 <- c(tmp4, seq(l))
                }
            }
        }
    }
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        Replicates=tmp4,
        MCC=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$Replicates <- as.character(df$Replicates)
    df$MCC[which(is.nan(df$MCC))] <- NA
    df$MCC[which(df$MCC < 0)] <- 0
    df
}

Plot_MCC <- function(method, E, outfile){
    df <- aggregateMCC(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=MCC))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="MCC", limits=c(0, 0.7))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

Plot_MCC_eachCCI <- function(method, E, outfile){
    df <- aggregateMCC_eachCCI(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=Replicates, fill=MCC))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~CCI_Complexity, scales="free_y", space="free_y")
    gg <- gg + scale_fill_viridis(name="MCC", limits=c(0, 0.7))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

aggregateMCC_merge <- function(E){
    BinMethods <- c(
        "labelpermutation",
        "labelpermutation2",
        "halpern",
        "cabelloaguilar",
        "previous_sctensor",
        "sctensor")
    df <- data.frame()
    for(i in seq_along(BinMethods)){
        tmp.df <- aggregateMCC_eachCCI(E, BinMethods[i])
        tmp.df <- cbind(tmp.df, BinMethods[i])
        df <- rbind(df, tmp.df)
    }
    colnames(df)[6] <- "BinMethods"
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$BinMethods <- factor(df$BinMethods, levels=BinMethods)
    df
}

Plot_MCC_Merge <- function(E, outfile){
    df <- aggregateMCC_merge(E=E)
    # Plot
    gg <- ggplot(df, aes(x=BinMethods, y=MCC, fill=BinMethods))
    gg <- gg + geom_boxplot()
    gg <- gg + facet_wrap(~CCI_Complexity)
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.title.x=element_blank())
    gg <- gg + theme(axis.text.x=element_blank())
    gg <- gg + theme(axis.title.y=element_text(size = 30))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(1, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=20, height=10)
}

# FPR
aggregateFPR <- function(E, Method){
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
                    Method, "/FPR/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                fpr <- unlist(fpr)
                target <- which(!is.nan(fpr))
                out = c(out, mean(fpr[target]))
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
        FPR=out)
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

aggregateFPR_eachCCI <- function(E, Method){
    out <- c()
    tmp1 <- c()
    tmp2 <- c()
    tmp3 <- c()
    tmp4 <- c()
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
                    Method, "/FPR/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                fpr <- unlist(fpr)
                target <- which(!is.nan(fpr))
                l <- length(target)
                out = c(out, fpr[target])
                tmp1 <- c(tmp1, rep(Celltypes[i], l))
                tmp2 <- c(tmp2, rep(CCIPatterns[j], l))
                tmp3 <- c(tmp3, rep(CCI_Complexity[k], l))
                if(l != 0){
                    tmp4 <- c(tmp4, seq(l))
                }
            }
        }
    }
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        Replicates=tmp4,
        FPR=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$Replicates <- as.character(df$Replicates)
    df
}

Plot_FPR <- function(method, E, outfile){
    df <- aggregateFPR(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=FPR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="FPR", limits=c(0, 0.09))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

Plot_FPR_eachCCI <- function(method, E, outfile){
    df <- aggregateFPR_eachCCI(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=Replicates, fill=FPR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~CCI_Complexity, scales="free_y", space="free_y")
    gg <- gg + scale_fill_viridis(name="FPR", limits=c(0, 0.09))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

aggregateFPR_merge <- function(E){
    BinMethods <- c(
        "labelpermutation",
        "labelpermutation2",
        "halpern",
        "cabelloaguilar",
        "previous_sctensor",
        "sctensor")
    df <- data.frame()
    for(i in seq_along(BinMethods)){
        tmp.df <- aggregateFPR_eachCCI(E, BinMethods[i])
        tmp.df <- cbind(tmp.df, BinMethods[i])
        df <- rbind(df, tmp.df)
    }
    colnames(df)[6] <- "BinMethods"
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$BinMethods <- factor(df$BinMethods, levels=BinMethods)
    df
}

Plot_FPR_Merge <- function(E, outfile){
    df <- aggregateFPR_merge(E=E)
    # Plot
    gg <- ggplot(df, aes(x=BinMethods, y=FPR, fill=BinMethods))
    gg <- gg + geom_boxplot()
    gg <- gg + facet_wrap(~CCI_Complexity)
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.title.x=element_blank())
    gg <- gg + theme(axis.text.x=element_blank())
    gg <- gg + theme(axis.title.y=element_text(size = 30))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(1, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=20, height=10)
}

# FNR
aggregateFNR <- function(E, Method){
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
                    Method, "/FNR/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                fnr <- unlist(fnr)
                target <- which(!is.nan(fnr))
                out = c(out, mean(fnr[target]))
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
        FNR=out)
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

aggregateFNR_eachCCI <- function(E, Method){
    out <- c()
    tmp1 <- c()
    tmp2 <- c()
    tmp3 <- c()
    tmp4 <- c()
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
                    Method, "/FNR/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                fnr <- unlist(fnr)
                target <- which(!is.nan(fnr))
                l <- length(target)
                out = c(out, fnr[target])
                tmp1 <- c(tmp1, rep(Celltypes[i], l))
                tmp2 <- c(tmp2, rep(CCIPatterns[j], l))
                tmp3 <- c(tmp3, rep(CCI_Complexity[k], l))
                if(l != 0){
                    tmp4 <- c(tmp4, seq(l))
                }
            }
        }
    }
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        Replicates=tmp4,
        FNR=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$Replicates <- as.character(df$Replicates)
    df
}

Plot_FNR <- function(method, E, outfile){
    df <- aggregateFNR(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=FNR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="FNR", limits=c(0, 1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

Plot_FNR_eachCCI <- function(method, E, outfile){
    df <- aggregateFNR_eachCCI(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=Replicates, fill=FNR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~CCI_Complexity, scales="free_y", space="free_y")
    gg <- gg + scale_fill_viridis(name="FNR", limits=c(0, 1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

aggregateFNR_merge <- function(E){
    BinMethods <- c(
        "labelpermutation",
        "labelpermutation2",
        "halpern",
        "cabelloaguilar",
        "previous_sctensor",
        "sctensor")
    df <- data.frame()
    for(i in seq_along(BinMethods)){
        tmp.df <- aggregateFNR_eachCCI(E, BinMethods[i])
        tmp.df <- cbind(tmp.df, BinMethods[i])
        df <- rbind(df, tmp.df)
    }
    colnames(df)[6] <- "BinMethods"
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$BinMethods <- factor(df$BinMethods, levels=BinMethods)
    df
}

Plot_FNR_Merge <- function(E, outfile){
    df <- aggregateFNR_merge(E=E)
    # Plot
    gg <- ggplot(df, aes(x=BinMethods, y=FNR, fill=BinMethods))
    gg <- gg + geom_boxplot()
    gg <- gg + facet_wrap(~CCI_Complexity)
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.title.x=element_blank())
    gg <- gg + theme(axis.text.x=element_blank())
    gg <- gg + theme(axis.title.y=element_text(size = 30))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(1, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=20, height=10)
}

# PR
aggregatePR <- function(E, Method){
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
                    Method, "/PR/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                pr <- unlist(pr)
                target <- which(!is.nan(pr))
                out = c(out, mean(pr[target]))
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
        PR=out)
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

aggregatePR_eachCCI <- function(E, Method){
    out <- c()
    tmp1 <- c()
    tmp2 <- c()
    tmp3 <- c()
    tmp4 <- c()
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
                    Method, "/PR/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                pr <- unlist(pr)
                target <- which(!is.nan(pr))
                l <- length(target)
                out = c(out, pr[target])
                tmp1 <- c(tmp1, rep(Celltypes[i], l))
                tmp2 <- c(tmp2, rep(CCIPatterns[j], l))
                tmp3 <- c(tmp3, rep(CCI_Complexity[k], l))
                if(l != 0){
                    tmp4 <- c(tmp4, seq(l))
                }
            }
        }
    }
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        Replicates=tmp4,
        PR=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$Replicates <- as.character(df$Replicates)
    df
}

Plot_PR <- function(method, E, outfile){
    df <- aggregatePR(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=PR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="PR", limits=c(0, 0.1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

Plot_PR_eachCCI <- function(method, E, outfile){
    df <- aggregatePR_eachCCI(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=Replicates, fill=PR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~CCI_Complexity, scales="free_y", space="free_y")
    gg <- gg + scale_fill_viridis(name="PR", limits=c(0, 0.1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

aggregatePR_merge <- function(E){
    BinMethods <- c(
        "labelpermutation",
        "labelpermutation2",
        "halpern",
        "cabelloaguilar",
        "previous_sctensor",
        "sctensor")
    df <- data.frame()
    for(i in seq_along(BinMethods)){
        tmp.df <- aggregatePR_eachCCI(E, BinMethods[i])
        tmp.df <- cbind(tmp.df, BinMethods[i])
        df <- rbind(df, tmp.df)
    }
    colnames(df)[6] <- "BinMethods"
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$BinMethods <- factor(df$BinMethods, levels=BinMethods)
    df
}

Plot_PR_Merge <- function(E, outfile){
    df <- aggregatePR_merge(E=E)
    # Plot
    gg <- ggplot(df, aes(x=BinMethods, y=PR, fill=BinMethods))
    gg <- gg + geom_boxplot()
    gg <- gg + facet_wrap(~CCI_Complexity)
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.title.x=element_blank())
    gg <- gg + theme(axis.text.x=element_blank())
    gg <- gg + theme(axis.title.y=element_text(size = 30))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(1, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=20, height=10)
}

# TR
TR <- function(trueCaH){
    tc <- lapply(trueCaH, function(t){
        length(which(t == 1)) / length(t)
    })
    sum(unlist(tc))
}

aggregateTR <- function(){
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
                inputfile <- paste0("data/groundtruth/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                out = c(out, TR(trueCaH))
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
        TR=out)
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

TR_eachCCI <- function(trueCaH){
    tc <- lapply(trueCaH, function(t){
        length(which(t == 1)) / length(t)
    })
    unlist(tc)
}

aggregateTR_eachCCI <- function(){
    out <- c()
    tmp1 <- c()
    tmp2 <- c()
    tmp3 <- c()
    tmp4 <- c()
    Celltypes=c("3","5","10","20","30")
    CCIPatterns=c("1","3","5")
    CCI_Complexity=c("OnetoOne", "ManytoMany")
    nCL <- length(Celltypes)
    nCP <- length(CCIPatterns)
    nCC <- length(CCI_Complexity)
    for(i in 1:nCL){
        for(j in 1:nCP){
            for(k in 1:nCC){
                inputfile <- paste0("data/groundtruth/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                out = c(out, TR_eachCCI(trueCaH))
                l <- length(trueCaH)
                tmp1 <- c(tmp1, rep(Celltypes[i], l))
                tmp2 <- c(tmp2, rep(CCIPatterns[j], l))
                tmp3 <- c(tmp3, rep(CCI_Complexity[k], l))
                if(l != 0){
                    tmp4 <- c(tmp4, seq(l))
                }
            }
        }
    }
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        Replicates=tmp4,
        TR=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$Replicates <- as.character(df$Replicates)
    df
}

Plot_TR <- function(outfile){
    df <- aggregateTR()
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=TR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="TR", limits=c(0, max(df$TR)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

Plot_TR_eachCCI <- function(outfile){
    df <- aggregateTR_eachCCI()
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=Replicates, fill=TR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~CCI_Complexity, scales="free_y", space="free_y")
    gg <- gg + scale_fill_viridis(name="TR", limits=c(0, max(df$TR)))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

# AUC
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

aggregateAUC_eachCCI <- function(E, Method){
    out <- c()
    tmp1 <- c()
    tmp2 <- c()
    tmp3 <- c()
    tmp4 <- c()
    Celltypes=c("3","5","10","20","30")
    CCIPatterns=c(1,3,5)
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
                l <- length(target)
                out = c(out, auc[target])
                tmp1 <- c(tmp1, rep(Celltypes[i], l))
                tmp2 <- c(tmp2, rep(CCIPatterns[j], l))
                tmp3 <- c(tmp3, rep(CCI_Complexity[k], l))
                if(l != 0){
                    tmp4 <- c(tmp4, seq(l))
                }
            }
        }
    }
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        Replicates=tmp4,
        AUC=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Replicates <- as.character(df$Replicates)
    df
}

Plot_AUC <- function(method, E, outfile){
    df <- aggregateAUC(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=AUC))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="AUC", limits=c(0, 1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

Plot_AUC_eachCCI <- function(method, E, outfile){
    df <- aggregateAUC_eachCCI(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=Replicates, fill=AUC))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~CCI_Complexity, scales="free_y", space="free_y")
    gg <- gg + scale_fill_viridis(name="AUC", limits=c(0, 1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

aggregateAUC_merge <- function(E){
    Methods <- c(
        "labelpermutation_tensor",
        "labelpermutation",
        "labelpermutation2_tensor",
        "labelpermutation2",
        "halpern_tensor",
        "halpern",
        "cabelloaguilar_tensor",
        "cabelloaguilar",
        "previous_sctensor",
        "sctensor")
    df <- data.frame()
    for(i in seq_along(Methods)){
        tmp.df <- aggregateAUC_eachCCI(E, Methods[i])
        tmp.df <- cbind(tmp.df, Methods[i])
        df <- rbind(df, tmp.df)
    }
    colnames(df)[6] <- "Methods"
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Methods <- factor(df$Methods, levels=Methods)
    df
}

Plot_AUC_Merge <- function(E, outfile){
    df <- aggregateAUC_merge(E=E)
    # Plot
    gg <- ggplot(df, aes(x=Methods, y=AUC, fill=Methods))
    gg <- gg + geom_boxplot()
    gg <- gg + facet_wrap(~CCI_Complexity)
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.title.x=element_blank())
    gg <- gg + theme(axis.text.x=element_blank())
    gg <- gg + theme(axis.title.y=element_text(size = 30))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(1, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=20, height=10)
}

# AUCPR
aggregateAUCPR <- function(E, Method){
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
                    Method, "/AUCPR/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                aucpr <- unlist(aucpr)
                target <- which(!is.nan(aucpr))
                out = c(out, mean(aucpr[target]))
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
        AUCPR=out)
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

aggregateAUCPR_eachCCI <- function(E, Method){
    out <- c()
    tmp1 <- c()
    tmp2 <- c()
    tmp3 <- c()
    tmp4 <- c()
    Celltypes=c("3","5","10","20","30")
    CCIPatterns=c(1,3,5)
    CCI_Complexity=c("OnetoOne", "ManytoMany")
    nCL <- length(Celltypes)
    nCP <- length(CCIPatterns)
    nCC <- length(CCI_Complexity)
    for(i in 1:nCL){
        for(j in 1:nCP){
            for(k in 1:nCC){
                inputfile <- paste0("output/",
                    E, "/",
                    Method, "/AUCPR/",
                    Celltypes[i], "Celltypes_",
                    CCIPatterns[j], "CCIPatterns_",
                    CCI_Complexity[k], ".RData")
                load(inputfile)
                aucpr <- unlist(aucpr)
                target <- which(!is.nan(aucpr))
                l <- length(target)
                out = c(out, aucpr[target])
                tmp1 <- c(tmp1, rep(Celltypes[i], l))
                tmp2 <- c(tmp2, rep(CCIPatterns[j], l))
                tmp3 <- c(tmp3, rep(CCI_Complexity[k], l))
                if(l != 0){
                    tmp4 <- c(tmp4, seq(l))
                }
            }
        }
    }
    df <- data.frame(
        Celltypes=tmp1,
        CCIPatterns=tmp2,
        CCI_Complexity=tmp3,
        Replicates=tmp4,
        AUCPR=out)
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Replicates <- as.character(df$Replicates)
    df
}

Plot_AUCPR <- function(method, E, outfile){
    df <- aggregateAUCPR(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=CCI_Complexity, fill=AUCPR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~., scales="free_y", space="free")
    gg <- gg + scale_fill_viridis(name="AUCPR", limits=c(0, 1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

Plot_AUCPR_eachCCI <- function(method, E, outfile){
    df <- aggregateAUCPR_eachCCI(E=E, Method=method)
    # Plot
    gg <- ggplot(df, aes(x=Celltypes, y=Replicates, fill=AUCPR))
    gg <- gg + geom_tile(color="white", size=0.1)
    gg <- gg + facet_grid(CCIPatterns~CCI_Complexity, scales="free_y", space="free_y")
    gg <- gg + scale_fill_viridis(name="AUCPR", limits=c(0, 1))
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + ylab("")
    gg <- gg + theme(axis.text.x=element_text(size = 30, hjust = 0), axis.title.x=element_blank())
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(3, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=13, height=10)
}

aggregateAUCPR_merge <- function(E){
    Methods <- c(
        "labelpermutation_tensor",
        "labelpermutation",
        "labelpermutation2_tensor",
        "labelpermutation2",
        "halpern_tensor",
        "halpern",
        "cabelloaguilar_tensor",
        "cabelloaguilar",
        "previous_sctensor",
        "sctensor")
    df <- data.frame()
    for(i in seq_along(Methods)){
        tmp.df <- aggregateAUCPR_eachCCI(E, Methods[i])
        tmp.df <- cbind(tmp.df, Methods[i])
        df <- rbind(df, tmp.df)
    }
    colnames(df)[6] <- "Methods"
    df$Celltypes <- factor(df$Celltypes,
        levels=unique(df$Celltypes))
    df$CCIPatterns <- factor(df$CCIPatterns,
        levels=rev(unique(df$CCIPatterns)))
    df$CCI_Complexity <- factor(df$CCI_Complexity,
        levels=c("OnetoOne", "ManytoMany"))
    df$Methods <- factor(df$Methods, levels=Methods)
    df
}

Plot_AUCPR_Merge <- function(E, outfile){
    df <- aggregateAUCPR_merge(E=E)
    # Plot
    gg <- ggplot(df, aes(x=Methods, y=AUCPR, fill=Methods))
    gg <- gg + geom_boxplot()
    gg <- gg + facet_wrap(~CCI_Complexity)
    gg <- gg + theme(strip.text = element_text(colour="blue3", face=2, size=40))
    gg <- gg + theme(axis.title.x=element_blank())
    gg <- gg + theme(axis.text.x=element_blank())
    gg <- gg + theme(axis.title.y=element_text(size = 30))
    gg <- gg + theme(axis.text.y=element_text(size = 30))
    gg <- gg + theme(legend.title = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.text = element_text(size = 30, hjust = 0))
    gg <- gg + theme(legend.key.size=unit(1, "cm"))
    gg <- gg + theme(legend.key.width=unit(1, "cm"))
    # Save
    ggsave(outfile, gg, dpi=120, width=20, height=10)
}


# ROC
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

Plot_PRC_AUCPR_MCC <- function(trueCaH, prc, aucpr, mcc, outfile){
    colvec <- c(brewer.pal(9, "Set1"),
        brewer.pal(8, "Set2"),
        brewer.pal(8, "Dark2"))
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
            cex=2, pch=16, type="b",
            xlab="Recall (Sensitivity)", ylab="Precision (PPV)",
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

Vec2Tensor_sctensor <- function(score, bin, truecah, ncelltypes){
    tmp <- unlist(lapply(score, function(s, t){
        AUC(as.vector(s@data), t)
    }, t=truecah))
    maxposition <- which(max(tmp) == tmp)[1]
    bin[[maxposition]]@data
}

Vec2Tensor <- function(bin, ncelltypes){
    counter <- 1
    npair <- length(bin) / ncelltypes^2
    out <- array(0, dim=c(ncelltypes, ncelltypes, npair))
    for(j in seq_len(npair)){
        tmp <- bin[counter:(counter+ncelltypes^2-1)]
        dim(tmp) <- c(ncelltypes, ncelltypes)
        out[,,j] <- tmp
        counter <- counter + ncelltypes^2
    }
    out
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

ROCCurve <- function(score, actual){
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
            list(auc=auc, roc=roc, f=f, aucpr=aucpr,
                prc=prc, mcc=mcc, fpr=fpr, fnr=fnr, pr=pr)
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

# 厳しいととってこない
# 1.5 => 1 => 0.8 => 0.5
.MAD <- function(x, thr=1.0){
    cluster <- abs(x - median(x))
    target <- which((x - median(x)) >= thr*median(cluster))
    cluster[target] <- "selected"
    cluster[setdiff(seq(cluster), target)] <- "not selected"
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

BIN_2 <- function(out, x){
    # L <- scTensor:::.HCLUST(out$ligand[x[1],])
    # R <- scTensor:::.HCLUST(out$receptor[x[2],])
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
