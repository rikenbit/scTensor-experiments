source("src/Functions.R")

# Parameter
sample <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]
# sample = "data/groundtruth/Mouse_Uterus.RData"
infile = paste0("data/groundtruth/", sample, ".RData")
load(infile)

l <- length(cif) - 1
png(file=outfile, width=mywidth(l), height=myheight(l))
mypanel(l)
lapply(2:(l+1),
	function(x){bipartiteGraph(cif, x)})
dev.off()

