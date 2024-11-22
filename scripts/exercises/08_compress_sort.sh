#!/usr/bin/env bash

# change directory to the alignments directory
cd ~/project/results/alignments

# compress the sam file
samtools view -bh SRR519926.sam > SRR519926.bam

# sort the bam file
samtools sort SRR519926.bam > SRR519926.sorted.bam

# index the sorted bam file
samtools index SRR519926.sorted.bam
