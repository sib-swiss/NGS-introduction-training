#!/usr/bin/env bash

WORKDIR=/config/workdir/project1
cd "$WORKDIR"/data/fastq

mkdir -p "$WORKDIR"/log

for SAMPLE in mother father son
do
    fastp \
    -o "$SAMPLE"_trimmed_R1.fastq.gz \
    -O "$SAMPLE"_trimmed_R2.fastq.gz \
    -i "$SAMPLE"_R1.fastq.gz \
    -I "$SAMPLE"_R2.fastq.gz \
    --detect_adapter_for_pe \
    --html "$SAMPLE".html \
    --json "$SAMPLE".json
done 