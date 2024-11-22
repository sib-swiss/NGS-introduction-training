#!/usr/bin/env bash

# change directory to the alignments directory
cd ~/project/results/alignments

# extract the unmapped reads with samtools view
# -b: output in BAM format
# -h: include the header in the output
# -f 0x4: only include reads that are unmapped
samtools view -bh -f 0x4 SRR519926.sorted.bam > SRR519926.sorted.unmapped.bam
