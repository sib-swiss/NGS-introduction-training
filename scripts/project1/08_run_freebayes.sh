#!/usr/bin/env bash

PROJDIR=/config/project/projects/project1
REFDIR="$PROJDIR"/data/reference/
ALIGNDIR="$PROJDIR"/alignments/
VCFDIR="$PROJDIR"/variants

mkdir -p "$VCFDIR"

samtools faidx "$REFDIR"/Homo_sapiens.GRCh38.dna.chromosome.20.fa

freebayes \
-f "$REFDIR"/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
-r chr20:10018000-10220000 \
"$ALIGNDIR"/merged.bam \
> "$VCFDIR"/variants.vcf
