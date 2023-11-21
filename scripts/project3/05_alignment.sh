#!/usr/bin/env bash

PROJDIR=/config/project/projects/project3/
cd "$PROJDIR"/data/reference

mkdir -p "$PROJDIR"/results/alignments

hisat2-build --threads 4 Mus_musculus.GRCm38.dna.chromosome.5.fa Mus_musculus.GRCm38.dna.chromosome.5.fa

for SAMPLE in Control1 Control2 Control3 Case1 Case2 Case3
do
    hisat2 \
    -x "$PROJDIR"/data/reference/Mus_musculus.GRCm38.dna.chromosome.5.fa \
    -1 "$PROJDIR"/results/trimmed/"$SAMPLE"_trimmed_R1.fastq.gz  \
    -2 "$PROJDIR"/results/trimmed/"$SAMPLE"_trimmed_R2.fastq.gz  \
    --threads 4 \
    | samtools sort \
    | samtools view -bh \
    > "$PROJDIR"/results/alignments/"$SAMPLE".bam

    samtools index "$PROJDIR"/results/alignments/"$SAMPLE".bam
done
