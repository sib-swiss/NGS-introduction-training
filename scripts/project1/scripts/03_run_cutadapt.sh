#!/usr/bin/env bash

WORKDIR=/config/workdir/project1
cd "$WORKDIR"/data/fastq

mkdir -p "$WORKDIR"/log

for SAMPLE in mother father son
do
    cutadapt \
    --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --quality-cutoff 10,10 \
    --minimum-length 25 \
    --output "$SAMPLE"_trimmed_R1.fastq.gz \
    --paired-output "$SAMPLE"_trimmed_R2.fastq.gz \
    "$SAMPLE"_R1.fastq.gz \
    "$SAMPLE"_R2.fastq.gz \
    > "$WORKDIR"/log/cutadapt_"$SAMPLE".log
done 