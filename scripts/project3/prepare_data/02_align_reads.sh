#!/usr/bin/env bash

#SBATCH --mail-user=geert.vangeest@bioinformatics.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-6
#SBATCH --job-name=align_reads
#SBATCH --partition=pall

PROJDIR=/data/users/gvangeest/courses/20230515_NGSQC

cd "$PROJDIR"/raw_data/
SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" sample_names.txt`

mkdir -p "$PROJDIR"/results/alignments/

apptainer exec \
/mnt/apps/centos7/hisat2_samtools_v4.0.0-beta.sif \
hisat2 \
-x /data/references/Mus_musculus/Ensembl/GRCm39/Sequence/hisat2/Mus_musculus.GRCm39.dna.primary_assembly.fa \
-1 "$SAMPLE"_S*_R1_001.fastq.gz \
-2 "$SAMPLE"_S*_R2_001.fastq.gz \
--rna-strandness RF \
-p "$SLURM_CPUS_PER_TASK" \
| samtools sort \
| samtools view -bh \
> "$PROJDIR"/results/alignments/"$SAMPLE".bam

samtools index "$PROJDIR"/results/alignments/"$SAMPLE".bam