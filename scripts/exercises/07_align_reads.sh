#!/usr/bin/env bash

# set the variables for the directory of trimmed reads, reference genome, and aligned reads
TRIMMED_DIR=~/project/results/trimmed
REFERENCE_DIR=~/project/ref_genome/
ALIGNED_DIR=~/project/results/alignments

# create the aligned directory which will contain the aligned reads
mkdir -p $ALIGNED_DIR

# align the trimmed reads to the reference genome using bowtie2
# -x: path to the bowtie2 index
# -1: path to the first trimmed read file
# -2: path to the second trimmed read file
bowtie2 \
-x $REFERENCE_DIR/ecoli-strK12-MG1655.fasta \
-1 $TRIMMED_DIR/trimmed_SRR519926_1.fastq \
-2 $TRIMMED_DIR/trimmed_SRR519926_2.fastq \
> $ALIGNED_DIR/SRR519926.sam
