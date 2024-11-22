#!/usr/bin/env bash

# create variables for directories
# TRIMMED_DIR is the directory where trimmed reads will be saved
# READS_DIR is the directory where the raw reads are saved
TRIMMED_DIR=~/project/results/trimmed
READS_DIR=~/project/reads

# create the trimmed directory
mkdir -p $TRIMMED_DIR

# change directory to the trimmed directory
cd $TRIMMED_DIR

# run fastp on the raw reads
# note that the backslash at the end of the line allows the command to continue on the next line
fastp \
-i $READS_DIR/SRR519926_1.fastq \
-I $READS_DIR/SRR519926_2.fastq \
-o $TRIMMED_DIR/trimmed_SRR519926_1.fastq \
-O $TRIMMED_DIR/trimmed_SRR519926_2.fastq \
--qualified_quality_phred 10 \
--length_required 25 \
--unqualified_percent_limit 80 \
--cut_front \
--cut_tail \
--detect_adapter_for_pe
