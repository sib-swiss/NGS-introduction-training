#!/usr/bin/env bash

PROJDIR=/config/project/projects/project1
REFDIR="$PROJDIR"/data/reference/
FASTQDIR="$PROJDIR"/data/fastq/
ALIGNDIR="$PROJDIR"/alignments/

mkdir -p "$ALIGNDIR"

cd "$REFDIR"

for SAMPLE in mother father son
do
    gatk MarkDuplicates \
    --INPUT "$ALIGNDIR"/"$SAMPLE".bam \
    --OUTPUT "$ALIGNDIR"/"$SAMPLE".md.bam \
    --METRICS_FILE "$ALIGNDIR"/"$SAMPLE".metrics.txt \
    2> "$PROJDIR"/log/markdup_"$SAMPLE".log

    samtools index "$ALIGNDIR"/"$SAMPLE".md.bam
done 
