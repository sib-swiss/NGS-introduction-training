#!/usr/bin/env bash

WORKDIR=/config/workdir/projects/project3/
cd "$WORKDIR"/data/

mkdir -p "$WORKDIR"/results/alignments

hisat2-build --threads 4 Mus_musculus.GRCm38.dna.chromosome.5.fa Mus_musculus.GRCm38.dna.chromosome.5.fa

for SAMPLE in Control1 Control2 Control3 Case1 Case2 Case3
do
    hisat2 \
    -x "$WORKDIR"/data/Mus_musculus.GRCm38.dna.chromosome.5.fa \
    -1 "$WORKDIR"/results/trimmed/"$SAMPLE"_trimmed_R1.fastq.gz  \
    -2 "$WORKDIR"/results/trimmed/"$SAMPLE"_trimmed_R2.fastq.gz  \
    --threads 4 \
    | samtools sort \
    | samtools view -bh \
    > "$WORKDIR"/results/alignments/"$SAMPLE".bam

    samtools index "$WORKDIR"/results/alignments/"$SAMPLE".bam
done
