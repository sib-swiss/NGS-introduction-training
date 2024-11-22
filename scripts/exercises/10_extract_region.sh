#!/usr/bin/env bash

# change directory to the alignments directory
cd ~/project/results/alignments

# extract a region from the sorted bam file
# U00096.3:2000000-2500000: region to extract
samtools view -bh \
SRR519926.sorted.bam \
U00096.3:2000000-2500000 \
> SRR519926.sorted.region.bam
