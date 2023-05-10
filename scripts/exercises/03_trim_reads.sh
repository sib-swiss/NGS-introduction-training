#!/usr/bin/env bash

TRIMMED_DIR=~/workdir/trimmed_data
READS_DIR=~/workdir/reads

mkdir -p $TRIMMED_DIR

fastp \


cutadapt \
--adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--quality-cutoff 10,10 \
--minimum-length 25 \
--output $TRIMMED_DIR/trimmed_SRR519926_1.fastq \
--paired-output $TRIMMED_DIR/trimmed_SRR519926_2.fastq \
$READS_DIR/SRR519926_1.fastq \
$READS_DIR/SRR519926_2.fastq
