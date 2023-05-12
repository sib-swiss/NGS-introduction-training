library(rtracklayer)
library(GenomicRanges)
gtf <- import.gff("ftp://ftp.ensembl.org/pub/release-102/gff3/mus_musculus/Mus_musculus.GRCm38.102.chromosome.5.gff3.gz")

gtf_sub <- subsetByOverlaps(gtf, as("5:100000000-120000000", "GRanges"),
                            type = "within")

export.gff2(gtf_sub, "results/Mus_musculus.GRCm38.102.chromosome.5.gtf")
