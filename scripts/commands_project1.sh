#!/bin/bash

### Download reads

mkdir -p project1/reads

cd project1/reads

wget -O sample_131_R1.fastq.gz \
https://genelab-data.ndc.nasa.gov/genelab/static/media/dataset/GLDS-251_rna-seq_13JUN2017HiSeq_Run_Sample_131_UMISS_Hoeksema_ACAGTG_L001_R1_001.fastq.gz?version=1

wget -O sample_131_R2.fastq.gz \
https://genelab-data.ndc.nasa.gov/genelab/static/media/dataset/GLDS-251_rna-seq_13JUN2017HiSeq_Run_Sample_131_UMISS_Hoeksema_ACAGTG_L001_R2_001.fastq.gz?version=1

### Download and index reference

cd ..

mkdir reference

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-48/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# for bowtie2
bowtie2-build --threads 4 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

# and hisat2
hisat2-build --threads 4 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

### Perform quality control

cd ../reads

# takes ~ 10 minutes:
fastqc --threads 2 *.fastq.gz

cd ..

TRIMMED_DIR=./trimmed_data
READS_DIR=./reads

mkdir $TRIMMED_DIR

# takes about 12 minutes on a single core:
cutadapt \
--adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--cores 4 \
--quality-cutoff 10,10 \
--minimum-length 25 \
--output $TRIMMED_DIR/trimmed_sample_131_R1.fastq.gz \
--paired-output  $TRIMMED_DIR/trimmed_sample_131_R2.fastq.gz \
$READS_DIR/sample_131_R1.fastq.gz \
$READS_DIR/sample_131_R2.fastq.gz

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
-x $REFERENCE_DIR/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
-1 $TRIMMED_DIR/trimmed_sample_131_R1.fastq.gz  \
-2 $TRIMMED_DIR/trimmed_sample_131_R2.fastq.gz  \
--threads 4 \
| samtools sort - \
| samtools view -bh - \
> $ALIGNED_DIR/sample_131.bt2.bam

samtools index $ALIGNED_DIR/sample_131.bt2.bam

# get only a subset of genes to view with IGV
samtools view -bh $ALIGNED_DIR/sample_131.bt2.bam 1:22145-42561 > $ALIGNED_DIR/sample_131.bt2.subset.bam
samtools index $ALIGNED_DIR/sample_131.bt2.subset.bam

# takes about 20 minutes:
hisat2 \
-x $REFERENCE_DIR/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
-1 $TRIMMED_DIR/trimmed_sample_131_R1.fastq.gz  \
-2 $TRIMMED_DIR/trimmed_sample_131_R2.fastq.gz  \
--threads 4 \
| samtools sort \
| samtools view -bh \
> $ALIGNED_DIR/sample_131.hs2.bam

samtools index $ALIGNED_DIR/sample_131.hs2.bam

# get only a subset of genes to view with IGV
samtools view -bh $ALIGNED_DIR/sample_131.hs2.bam 1:22145-42561 > $ALIGNED_DIR/sample_131.hs2.subset.bam
samtools index $ALIGNED_DIR/sample_131.hs2.subset.bam

### count features

cd ./reference
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GTF_genes_transposons.Mar202021.gtf.gz
gunzip Araport11_GTF_genes_transposons.Mar202021.gtf.gz
sed 's/^Chr//g' Araport11_GTF_genes_transposons.Mar202021.gtf > Araport11_GTF_genes_transposons.Mar202021.noChr.gtf

cd ..

mkdir ./counts

featureCounts \
-p \
-T 2 \
-a $REFERENCE_DIR/Araport11_GTF_genes_transposons.Mar202021.noChr.gtf \
-o ./counts/counts.txt \
$ALIGNED_DIR/sample_131.bt2.bam $ALIGNED_DIR/sample_131.hs2.bam
