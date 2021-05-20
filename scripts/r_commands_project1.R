cts <- read.delim("project_work/project1/counts/counts.txt", comment.char = "#", row.names = 1)

cts <- cts[,c(6,7)]
colnames(cts) <- gsub('..alignments.', '', colnames(cts))

logcpm <- apply(cts, 2, function(x){
  cpm <- x/sum(x)*1e6
  out <- log2(cpm + 0.1)
  return(out)
})
logcpm <- as.data.frame(logcpm)

bt2 <- logcpm[,"sample_131.bt2.bam"]
hs2 <- logcpm[,"sample_131.hs2.bam"]


plot(bt2,hs2, pch = 19, cex = 0.3, xlab = "bowtie2", ylab = "hisat2")
abline(0,1, col = 'red')

diff <- abs(bt2 - hs2)
names(diff) <- rownames(cts)
top10 <- diff[order(diff, decreasing = TRUE)][1:10]
View(cts[names(top10),])


# check out e.g. AT1G59359 (low with hisat2): most reads have their mate in another annotation
# check out e.g. AT1G65970 (high with hisat2): 