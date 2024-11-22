#!/usr/bin/env bash

# change directory to the reference directory
cd ~/project/ref_genome

# build the bowtie2 index
# note that bowtie2-build takes two positional arguments
# the first is the reference genome file
# the second is the basename of the index files 
# these are typically the same by convention
bowtie2-build ecoli-strK12-MG1655.fasta ecoli-strK12-MG1655.fasta
