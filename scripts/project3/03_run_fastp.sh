#!/usr/bin/env bash

PROJDIR=/config/project/projects/project3/
cd "$PROJDIR"/data/fastq

mkdir -p "$PROJDIR"/results/trimmed
mkdir -p "$PROJDIR"/log

for SAMPLE in Control1 Control2 Control3 Case1 Case2 Case3
do
    fastp \
    -o "$PROJDIR"/results/trimmed/"$SAMPLE"_trimmed_R1.fastq.gz \
    -O "$PROJDIR"/results/trimmed/"$SAMPLE"_trimmed_R2.fastq.gz \
    -i "$SAMPLE"_R1.fastq.gz \
    -I "$SAMPLE"_R2.fastq.gz \
    --detect_adapter_for_pe \
    --html "$PROJDIR"/results/trimmed/"$SAMPLE".html \
    --json "$PROJDIR"/results/trimmed/"$SAMPLE".json
done 
