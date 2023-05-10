#!/usr/bin/env bash

WORKDIR=/config/workdir/project1
REFDIR="$WORKDIR"/data/reference/
FASTQDIR="$WORKDIR"/data/fastq/
ALIGNDIR="$WORKDIR"/data/alignments/

mkdir -p "$ALIGNDIR"

cd "$REFDIR"

bowtie2-build Homo_sapiens.GRCh38.dna.chromosome.20.fa \
Homo_sapiens.GRCh38.dna.chromosome.20.fa

for SAMPLE in mother father son
do
    bowtie2 \
    -x $REFDIR/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    -1 $FASTQDIR/"$SAMPLE"_trimmed_R1.fastq.gz \
    -2 $FASTQDIR/"$SAMPLE"_trimmed_R2.fastq.gz \
    --rg-id "$SAMPLE" \
    --rg SM:"$SAMPLE" \
    2> "$WORKDIR"/log/bowtie2_"$SAMPLE".log \
    | samtools sort - \
    | samtools view -bh - \
    > "$ALIGNDIR"/"$SAMPLE".bam

    samtools index "$ALIGNDIR"/"$SAMPLE".bam
done 
