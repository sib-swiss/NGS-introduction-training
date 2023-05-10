#!/usr/bin/env bash

WORKDIR=/config/workdir/project1
REFDIR="$WORKDIR"/data/reference/
ALIGNDIR="$WORKDIR"/data/alignments/
VCFDIR="$WORKDIR"/data/variants

freebayes \
-f "$REFDIR"/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
-r chr20:10018000-10220000 \
"$ALIGNDIR"/merged.bam \
2> "$WORKDIR"/log/freebayes.log \
> "$VCFDIR"/variants.vcf
