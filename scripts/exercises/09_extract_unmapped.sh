#!/usr/bin/env bash

cd ~/project/results/alignments

samtools view -bh -f 0x4 SRR519926.sorted.bam > SRR519926.sorted.unmapped.bam
