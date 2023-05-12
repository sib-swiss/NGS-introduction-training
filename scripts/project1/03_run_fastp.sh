#!/usr/bin/env bash

WORKDIR=/config/workdir/projects/project1
cd "$WORKDIR"/data/fastq

mkdir -p "$WORKDIR"/results/trimmed
mkdir -p "$WORKDIR"/log

for SAMPLE in mother father son
do
    fastp \
    -o "$WORKDIR"/results/trimmed/"$SAMPLE"_trimmed_R1.fastq.gz \
    -O "$WORKDIR"/results/trimmed/"$SAMPLE"_trimmed_R2.fastq.gz \
    -i "$SAMPLE"_R1.fastq.gz \
    -I "$SAMPLE"_R2.fastq.gz \
    --detect_adapter_for_pe \
    --html "$WORKDIR"/results/trimmed/"$SAMPLE".html \
    --json "$WORKDIR"/results/trimmed/"$SAMPLE".json
done 
