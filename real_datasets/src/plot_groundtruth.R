source("src/Functions.R")

# Parameter
sample <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]
# sample = "Human_NicotinehESCs_Nicotine"
infile = paste0("data/groundtruth/", sample, ".RData")
load(infile)

l <- length(cif) - 1
png(file=outfile, width=mywidth(l), height=myheight(l))
mypanel(l)
lapply(seq(l), function(x){bipartiteGraph(cif, x+1)})
dev.off()


