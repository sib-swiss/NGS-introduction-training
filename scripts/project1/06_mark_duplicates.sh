#!/usr/bin/env bash

WORKDIR=/config/workdir/projects/project1
REFDIR="$WORKDIR"/data/reference/
FASTQDIR="$WORKDIR"/data/fastq/
ALIGNDIR="$WORKDIR"/alignments/

mkdir -p "$ALIGNDIR"

cd "$REFDIR"

for SAMPLE in mother father son
do
    gatk MarkDuplicates \
    --INPUT "$ALIGNDIR"/"$SAMPLE".bam \
    --OUTPUT "$ALIGNDIR"/"$SAMPLE".md.bam \
    --METRICS_FILE "$ALIGNDIR"/"$SAMPLE".metrics.txt \
    2> "$WORKDIR"/log/markdup_"$SAMPLE".log

    samtools index "$ALIGNDIR"/"$SAMPLE".md.bam
done 
