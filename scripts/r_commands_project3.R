cts <- read.delim("project_work/project3/counts/counts.txt", comment.char = "#", row.names = 1)

cts <- cts[,c(6,7)]
colnames(cts) <- gsub('..alignments.', '', colnames(cts))

logcpm <- apply(cts, 2, function(x){
  cpm <- x/sum(x)*1e6
  out <- log2(cpm + 0.1)
  return(out)
})
logcpm <- as.data.frame(logcpm)

bt2 <- logcpm$SRR7822040.chr5.bt2.bam
hs2 <- logcpm$SRR7822040.chr5.hs2.bam


plot(bt2,hs2, pch = 19, cex = 0.3, xlab = "bowtie2", ylab = "hisat2")
abline(0,1, col = 'red')

diff <- abs(bt2 - hs2)
names(diff) <- rownames(cts)
top10 <- diff[order(diff, decreasing = TRUE)][1:10]
View(cts[names(top10),])

# check out ENSMUSG00000106930
# look at an aligned read with bt2 (and not aligned with hs2) and check where it is in the bam file of hs2
# e.g. samtools view SRR7822040.chr5.hs2.bam | grep SRR7822040.17131809 | cut -f 1-5

