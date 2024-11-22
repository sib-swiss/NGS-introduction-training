#!/usr/bin/env bash

# create a variable for the reference directory
REFERENCE_DIR=~/project/ref_genome/

# create the reference directory
mkdir $REFERENCE_DIR
# change directory to the reference directory
cd $REFERENCE_DIR

# download the E. coli K-12 MG1655 reference genome
esearch -db nuccore -query 'U00096' \
| efetch -format fasta > ecoli-strK12-MG1655.fasta