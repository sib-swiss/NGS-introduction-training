library(DESeq2)

# function to read in and simplify column names of counts matrix
read_count_matrix <- function(counts_file) {
  cts <- read.delim(counts_file, comment.char = "#")
  
  # get columns that end with .bam
  sample_cols <- grepl("bam$", colnames(cts))

  # split by dot (/ of paths are replaced by dots)
  # and take the second last string (string before .bam)
  colnames(cts)[sample_cols] <- sapply(colnames(cts)[sample_cols],
                                       function(x) {
                                         spl <- strsplit(x, "\\.")[[1]]
                                         return(spl[length(spl) - 1])
                                       })
  
  # generate a count matrix with samples in columns and genes in rows
  count_matrix <- as.matrix(cts[,sample_cols])
  rownames(count_matrix) <- cts$Geneid
  return(count_matrix)
}

# apply function
count_matrix <- read_count_matrix("projects/project3/results/counts/counts.txt")

# create metadata object that describes the samples
col_data <- data.frame(condition = gsub("[0-9]", "", colnames(count_matrix)))
rownames(col_data) <- colnames(count_matrix)
View(col_data)

# run standard/default function for differential gene expression analysis
# according to https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd)
