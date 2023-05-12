#!/usr/bin/env bash

WORKDIR=/config/workdir/projects/project3/
cd "$WORKDIR"/data/

mkdir -p "$WORKDIR"/results/counts

featureCounts \
-p \
-T 4 \
-g Name \
-a Mus_musculus.GRCm38.102.chromosome.5.gtf \
-o "$WORKDIR"/results/counts/counts.txt \
"$WORKDIR"/results/alignments/*.bam
