#!/usr/bin/env bash

cd ~/workspace/projects/project2

mkdir -p alignments

for sample in EV_2 Cell_2; do
    minimap2 \
    -a \
    -x splice \
    -t 4 \
    references/Homo_sapiens.GRCh38.dna.primary_assembly.chr5.chr6.chrX.fa \
    reads/"$sample".fastq.gz \
    | samtools sort \
    | samtools view -bh > alignments/"$sample".bam

    ## indexing for IGV
    samtools index alignments/"$sample".bam
done