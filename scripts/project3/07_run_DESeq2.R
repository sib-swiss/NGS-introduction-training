library(DESeq2)

read_count_matrix <- function(counts_file) {
  cts <- read.delim(counts_file, comment.char = "#")
  
  sample_cols <- grepl("bam$", colnames(cts))
  colnames(cts)[sample_cols] <- sapply(colnames(cts)[sample_cols],
                                       function(x) {
                                         spl <- strsplit(x, "\\.")[[1]]
                                         return(spl[length(spl) - 1])
                                       })
  
  count_matrix <- as.matrix(cts[,sample_cols])
  rownames(count_matrix) <- cts$Geneid
  return(count_matrix)
}

count_matrix <- read_count_matrix("projects/project3/results/counts/counts.txt")

col_data <- data.frame(condition = gsub("[0-9]", "", colnames(count_matrix)))
rownames(col_data) <- colnames(count_matrix)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd)
