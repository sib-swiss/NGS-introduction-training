#!/usr/bin/env bash

TRIMMED_DIR=~/workdir/trimmed_data
REFERENCE_DIR=~/workdir/ref_genome
ALIGNED_DIR=~/workdir/alignment_output

bowtie2 \
-x $REFERENCE_DIR/ecoli-strK12-MG1655.fasta \
-1 $TRIMMED_DIR/paired_trimmed_SRR519926_1.fastq \
-2 $TRIMMED_DIR/paired_trimmed_SRR519926_2.fastq \
| samtools sort - \
| samtools view -bh - \
> $ALIGNED_DIR/SRR519926.sorted.mapped.frompipe.bam
