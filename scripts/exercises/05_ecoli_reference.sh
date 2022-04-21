#!/usr/bin/env bash

REFERENCE_DIR=~/workdir/ref_genome/

mkdir $REFERENCE_DIR
cd $REFERENCE_DIR

esearch -db nuccore -query 'U00096' \
| efetch -format fasta > ecoli-strK12-MG1655.fasta