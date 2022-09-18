source("src/Functions.R")

# Parameter
# infile1 = "data/Human_Germline_Female/Human_Germline_Female.RData"
# infile2 = "output/sctensor/Human_Germline_Female.RData"

infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Data loading
load(infile1)
load(infile2)

# Preprocessing
sce <- SingleCellExperiment(assays=list(counts = input))
DATA <- cbind(LR, NA, "DLRP")
colnames(DATA)[3:4] <- c("SOURCEID", "SOURCEDB")
if(length(grep("Human_", infile1)) != 0){
    METADATA <- data.frame(NAME="TAXID", VALUE="9606")
}else{
    METADATA <- data.frame(NAME="TAXID", VALUE="10090")
}

# SQLite setting
dbfile <- tempfile()
mydb <- dbConnect(SQLite(), dbfile)
dbWriteTable(mydb, "DATA", DATA)
dbWriteTable(mydb, "METADATA", METADATA)

# SingleCellExperiment setting
out.umap <- umap(t(input))
reducedDims(sce) <- SimpleList(UMAP=out.umap)

# Metadata setting
if(length(grep("Human_", infile1)) != 0){
    ahid <- "AH97772"
}else{
    ahid <- "AH97793"
}
metadata(sce) <- list(lrbase=dbfile,
    ahid=ahid,
    lr.evidence="known",
    color=celltypes,
    label=names(celltypes),
    algorithm="ntd2",
    sctensor=out,
    ranks=c(nrow(out$ligand), nrow(out$receptor)),
    datasize=c(ncol(out$ligand), ncol(out$receptor),
        dim(out$lrpair)[3]),
    recerror=out$recerror,
    relchange=out$relchange)

# HTML Report
setAnnotationHubOption("CACHE", getwd())
out.dir <- gsub("/index.html", "", outfile)
cellCellReport(sce,
    reducedDimNames="UMAP",
    out.dir=out.dir,
    html.open=FALSE,
    title="The result of scTensor",
    author="Koki Tsuyuzaki",
    assayNames="counts",
    thr=100,
    top="full",
    p=0.05,
    upper=10,
    goenrich=TRUE,
    meshenrich=TRUE,
    reactomeenrich=TRUE,
    doenrich=TRUE,
    ncgenrich=TRUE,
    dgnenrich=TRUE,
    nbins=40)
