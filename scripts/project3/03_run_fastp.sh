#!/usr/bin/env bash

WORKDIR=/config/workdir/projects/project3/
cd "$WORKDIR"/data/fastq

mkdir -p "$WORKDIR"/results/trimmed
mkdir -p "$WORKDIR"/log

for SAMPLE in Control1 Control2 Control3 Case1 Case2 Case3
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
