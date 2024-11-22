#!/usr/bin/env bash

# set the directories for the trimmed reads, reference genome, and aligned reads
TRIMMED_DIR=~/project/results/trimmed
REFERENCE_DIR=~/project/ref_genome
ALIGNED_DIR=~/project/results/alignments

# align the trimmed reads to the reference genome using bowtie2, sort the output, and save as a bam file
# all in one pipe
# -x: path to the bowtie2 index
# -1: path to the first trimmed read file
# -2: path to the second trimmed read file
# 2>: redirect stderr to a log file
# pipe the output of bowtie2 to samtools sort
# pipe the output of samtools sort to samtools view
# -: read from stdin
# >: redirect the output to a bam file
bowtie2 \
-x $REFERENCE_DIR/ecoli-strK12-MG1655.fasta \
-1 $TRIMMED_DIR/trimmed_SRR519926_1.fastq \
-2 $TRIMMED_DIR/trimmed_SRR519926_2.fastq \
2> $ALIGNED_DIR/bowtie2_SRR519926.log \
| samtools sort - \
| samtools view -bh - \
> $ALIGNED_DIR/SRR519926.sorted.mapped.frompipe.bam
