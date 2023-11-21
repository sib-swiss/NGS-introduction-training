#!/usr/bin/env bash

PROJDIR=/config/project/projects/project3/
cd "$PROJDIR"/data/reference

mkdir -p "$PROJDIR"/results/counts

featureCounts \
-p \
-T 4 \
-a Mus_musculus.GRCm38.102.chromosome.5.gtf \
-o "$PROJDIR"/results/counts/counts2.txt \
"$PROJDIR"/results/alignments/*.bam
