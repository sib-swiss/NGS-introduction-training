#!/usr/bin/env bash

#SBATCH --mail-user=geert.vangeest@bioinformatics.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-6
#SBATCH --job-name=align_reads
#SBATCH --partition=pall

PROJDIR=/data/users/gvangeest/courses/20230515_NGSQC

mkdir -p "$PROJDIR"/results/reads/

cd "$PROJDIR"/raw_data/
SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" sample_names.txt`

samtools view -bh "$PROJDIR"/results/alignments/"$SAMPLE".bam \
5:100000000-120000000 \
> "$PROJDIR"/results/alignments/"$SAMPLE"_chr5.bam

samtools collate \
-O \
-@ "$SLURM_CPUS_PER_TASK" \
"$PROJDIR"/results/alignments/"$SAMPLE"_chr5.bam \
"$PROJDIR"/results/alignments/tmp."$SAMPLE" \
| samtools fastq \
-@ "$SLURM_CPUS_PER_TASK" \
-1 "$PROJDIR"/results/reads/"$SAMPLE"_R1.fastq.gz \
-2 "$PROJDIR"/results/reads/"$SAMPLE"_R2.fastq.gz \
-s "$PROJDIR"/results/reads/"$SAMPLE"_S.fastq.gz \
-0 "$PROJDIR"/results/reads/"$SAMPLE"_U.fastq.gz

cat filename_conversions.txt | while read OLD NEW
do
    mv "$OLD"_R1.fastq.gz "$NEW"_R1.fastq.gz
    mv "$OLD"_R2.fastq.gz "$NEW"_R2.fastq.gz
done