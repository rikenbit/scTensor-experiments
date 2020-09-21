# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1] 
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
infile4 <- commandArgs(trailingOnly=TRUE)[4]
infile5 <- commandArgs(trailingOnly=TRUE)[5]
outfile <- commandArgs(trailingOnly=TRUE)[6]

# Data Loading
rss <- c()
load(infile1)
rss <- c(rss, rev(out$recerror)[1])
load(infile2)
rss <- c(rss, rev(out$recerror)[1])
load(infile3)
rss <- c(rss, rev(out$recerror)[1])
load(infile4)
rss <- c(rss, rev(out$recerror)[1])
load(infile5)
rss <- c(rss, rev(out$recerror)[1])
target <- which(min(rss) == rss)

# Copy
bestfitfile <- eval(parse(text=paste0("infile", target)))
cmd <- paste("system('cp -rf", bestfitfile, outfile, "')")
eval(parse(text=cmd))
