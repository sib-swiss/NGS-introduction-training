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
/mnt/apps/centos7/fastqc_0.11.9--hdfd78af_1.sif \
fastqc "$SAMPLE"_R?.fastq.gz
