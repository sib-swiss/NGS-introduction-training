#!/usr/bin/env bash

PROJDIR=/config/project/projects/project1
REFDIR="$PROJDIR"/data/reference/
FASTQDIR="$PROJDIR"/results/trimmed/
ALIGNDIR="$PROJDIR"/alignments/

mkdir -p "$PROJDIR"/log/

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
    2> "$PROJDIR"/log/bowtie2_"$SAMPLE".log \
    | samtools sort - \
    | samtools view -bh - \
    > "$ALIGNDIR"/"$SAMPLE".bam

    samtools index "$ALIGNDIR"/"$SAMPLE".bam
done 
