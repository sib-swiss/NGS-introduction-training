#!/usr/bin/env bash

WORKDIR=/config/workdir/projects/project1
REFDIR="$WORKDIR"/data/reference/
ALIGNDIR="$WORKDIR"/alignments/
VCFDIR="$WORKDIR"/variants

mkdir -p "$VCFDIR"

samtools faidx "$REFDIR"/Homo_sapiens.GRCh38.dna.chromosome.20.fa

freebayes \
-f "$REFDIR"/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
-r chr20:10018000-10220000 \
"$ALIGNDIR"/merged.bam \
> "$VCFDIR"/variants.vcf
