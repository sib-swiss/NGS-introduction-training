
PROJDIR=/data/users/gvangeest/courses/20230515_NGSQC

mkdir -p "$PROJDIR"/raw_data
cd "$PROJDIR"/raw_data

cat "$PROJDIR"/NGS-introduction-training/scripts/project3/prepare_data/s3_uris_project3.txt \
| while read DATE TIME SIZE FILE
do
    echo downloading $FILE
    wget https://sra-data-delivery-sib.s3.amazonaws.com/"$FILE" .
done
# s3://sra-pub-src-3/SRR7821918/TEI696A32_S1_L001_R1_001.fastq.gz