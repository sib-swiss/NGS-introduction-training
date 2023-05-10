#!/usr/bin/env bash

#SBATCH --mail-user=geert.vangeest@bioinformatics.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-6
#SBATCH --job-name=align_reads
#SBATCH --partition=pall

PROJDIR=/data/users/gvangeest/courses/20230515_NGSQC

cd "$PROJDIR"/raw_data/
SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" sample_names.txt`

cd "$PROJDIR"/results/reads/

singularity exec \
/mnt/apps/centos7/fastp_0.23.2--h5f740d0_3.sif \
fastp \
-i "$SAMPLE"_R1.fastq.gz \
-I "$SAMPLE"_R2.fastq.gz \
-o "$SAMPLE"_R1_trimmed.fastq.gz \
-O "$SAMPLE"_R2_trimmed.fastq.gz \
-h fastp_report_"$SAMPLE".html \
-j fastp_report_"$SAMPLE".json


