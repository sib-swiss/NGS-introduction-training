library(rtracklayer)
library(GenomicRanges)
gtf <- import.gff("ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz")

gtf_sub <- subsetByOverlaps(gtf, as("5:100000000-120000000", "GRanges"))

export.gff3(gtf_sub, "projects/project3/data/reference/Mus_musculus.GRCm38.102.chromosome.5.gtf")
