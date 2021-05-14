#!/bin/bash

### Download reads

mkdir -p project3/reads

cd project3/reads

wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project3/SRR7822040.chr5_R1.fastq.gz
wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project3/SRR7822040.chr5_R2.fastq.gz


### Download and index reference

cd ..

mkdir reference

wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project3/Mus_musculus.GRCm38.dna.primary_assembly.chr5.fa.gz

gunzip Mus_musculus.GRCm38.dna.primary_assembly.chr5.fa.gz

cd reference

# for bowtie2
bowtie2-build --threads 4 Mus_musculus.GRCm38.dna.primary_assembly.chr5.fa Mus_musculus.GRCm38.dna.primary_assembly.chr5.fa

# and hisat2
hisat2-build --threads 4 Mus_musculus.GRCm38.dna.primary_assembly.chr5.fa Mus_musculus.GRCm38.dna.primary_assembly.chr5.fa

### Perform quality control

cd ../reads

# takes ~ 10 minutes:
fastqc --threads 2 *.fastq.gz

cd ..

TRIMMED_DIR=./trimmed_data
READS_DIR=./reads

mkdir $TRIMMED_DIR

cutadapt \
--adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--cores 4 \
--quality-cutoff 10,10 \
--minimum-length 25 \
--output $TRIMMED_DIR/paired_trimmed_SRR7822040.chr5_R1.fastq.gz \
--paired-output  $TRIMMED_DIR/paired_trimmed_SRR7822040.chr5_R2.fastq.gz \
$READS_DIR/SRR7822040.chr5_R1.fastq.gz \
$READS_DIR/SRR7822040.chr5_R2.fastq.gz

cd $TRIMMED_DIR

fastqc --threads 2 *.fastq.gz

#https://sequencing.qcfail.com/articles/positional-sequence-bias-in-random-primed-libraries/

### Align reads

cd ..

READS_DIR=./reads
REFERENCE_DIR=./reference
ALIGNED_DIR=./alignments

mkdir -p $ALIGNED_DIR

# takes about 40 minutes
bowtie2 \
-x $REFERENCE_DIR/Mus_musculus.GRCm38.dna.primary_assembly.chr5.fa \
-1 $TRIMMED_DIR/paired_trimmed_SRR7822040.chr5_R1.fastq.gz  \
-2 $TRIMMED_DIR/paired_trimmed_SRR7822040.chr5_R2.fastq.gz  \
--threads 4 \
| samtools sort - \
| samtools view -bh - \
> $ALIGNED_DIR/SRR7822040.chr5.bt2.bam

samtools index $ALIGNED_DIR/SRR7822040.chr5.bt2.bam

# get only a subset of genes to view with IGV
samtools view -bh $ALIGNED_DIR/SRR7822040.chr5.bt2.bam 5:32592000-32999545 > $ALIGNED_DIR/SRR7822040.chr5.bt2.subset.bam
samtools index $ALIGNED_DIR/SRR7822040.chr5.bt2.subset.bam

# takes about 20 minutes:
hisat2 \
-x $REFERENCE_DIR/Mus_musculus.GRCm38.dna.primary_assembly.chr5.fa \
-1 $TRIMMED_DIR/paired_trimmed_SRR7822040.chr5_R1.fastq.gz  \
-2 $TRIMMED_DIR/paired_trimmed_SRR7822040.chr5_R2.fastq.gz  \
--threads 4 \
| samtools sort \
| samtools view -bh \
> $ALIGNED_DIR/SRR7822040.chr5.hs2.bam

samtools index $ALIGNED_DIR/SRR7822040.chr5.hs2.bam

# get only a subset of genes to view with IGV
samtools view -bh $ALIGNED_DIR/SRR7822040.chr5.hs2.bam 5:32592000-32999545 > $ALIGNED_DIR/SRR7822040.chr5.hs2.subset.bam
samtools index $ALIGNED_DIR/SRR7822040.chr5.hs2.subset.bam

### count features

wget http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz
gunzip mus_musculus/Mus_musculus.GRCm39.104.gtf.gz

cd ..

mkdir ./counts

featureCounts \
-p \
-T 2 \
-a $REFERENCE_DIR/Mus_musculus.GRCm38.102.chr5.gtf \
-o ./counts/counts.txt \
$ALIGNED_DIR/SRR7822040.chr5.bt2.bam $ALIGNED_DIR/SRR7822040.chr5.hs2.bam
