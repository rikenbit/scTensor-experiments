source("src/Functions.R")

# Parameter
sample <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]
# sample = "data/groundtruth/30Celltypes_5CCIPatterns_ManytoMany.RData"
infile = paste0("data/groundtruth/", sample, ".RData")
load(infile)

l <- length(cif) - 1
png(file=outfile, width=300*l, height=300)
layout(t(seq(l)))
lapply(2:(l+1),
	function(x){bipartiteGraph(cif, x)})
dev.off()
