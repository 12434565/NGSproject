dir="/wellstein/user_homes/il1380/project"
ref=${dir}/data/ref
tar -xzvf selected_dataset_files.tgz

# step1 qc
mkdir /wellstein/user_homes/il1380/project/1qc
cd /wellstein/user_homes/il1380/project/1qc
fastqc ../data/*.fastq.gz -o ${dir}/1qc -t 8
mv *_fastqc* ${dir}/1qc
multiqc ${dir}/1qc -o ${dir}/1qc/multiqc


# step2 Quality and adapter trimming and filtering
mkdir ${dir}/2trim
n=1
for file in ${dir}/data/*.fastq.gz; do
  sample=$(basename $file .fastq.gz)
  fastp -i $file \
        -o ${dir}/2trim/${sample}_trimmed.fastq.gz \
        -h ${dir}/2trim/${sample}_fastp.html \
        -j ${dir}/2trim/${sample}_fastp.json \
        -q 20 \
        -l 20 \
        --trim_front1 16 \
        -w 8
    echo "\nProcessed $n: $file\n"
    n=$((n + 1))
done
mkdir ${dir}/2trim/multiqc3
# multiqc ${dir}/2trim -o ${dir}/2trim/multiqc3
# mkdir ${dir}/2trim/fastqc
fastqc ${dir}/2trim/*.fastq.gz -o ${dir}/2trim/fastqc -t 8
multiqc ${dir}/2trim/fastqc -o ${dir}/2trim/multiqc3
echo "done" && python3 ${dir}/code/a.py

# step3 Alignment with GRCh38
# 3.1 download reference genome and annotation
cd /wellstein/user_homes/il1380/project/data/ref
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.gtf.gz
gunzip GCF_000001405.26_GRCh38_genomic.fna.gz
gunzip GCF_000001405.26_GRCh38_genomic.gtf.gz
hisat2-build GCF_000001405.26_GRCh38_genomic.fna hg38
echo "done" && python3 ${dir}/code/a.py
# screen -S hisat
screen -R hisat

# step3 align
mkdir -p ${dir}/3align

for file in ${dir}/2trim/*_trimmed.fastq.gz; do
  sample=$(basename $file _trimmed.fastq.gz)

  hisat2 -p 8 \
    -x ${ref}/hg38 \
    -U ${file} \
  | samtools view -bS - \
  | samtools sort -o ${dir}/3align/${sample}.bam
  samtools index ${dir}/3align/${sample}.bam

  echo "Aligned: $sample" 
done
echo "done" && python3 ${dir}/code/a.py

echo -e "Sample\tTotal\tMapped\tMappingRate" > alignment_summary.txt


# 获得统计值
# echo "" > flagstat_all.txt
# for bam in ${dir}/3align/*.bam
# do
#   echo "===== $(basename $bam) =====" >> flagstat_all.txt
#   samtools flagstat $bam >> flagstat_all.txt
#   echo "" >> flagstat_all.txt
# done
# 先清空/创建总log文件
all_log=${dir}/3align/hisat2_all.log
> $all_log

for file in ${dir}/2trim/*_trimmed.fastq.gz; do
  sample=$(basename $file _trimmed.fastq.gz)
  echo "===== $sample =====" >> $all_log
  hisat2 -p 8 \
    -x ${ref}/hg38 \
    -U ${file} \
    -S /dev/null \
    2>> $all_log
  echo "" >> $all_log
  echo "Aligned (log only): $sample"
done
echo "done" && python3 ${dir}/code/a.py

echo "done"
# step4 count
mkdir -p ${dir}/4counts

featureCounts -T 8 \
  -a ${ref}/GCF_000001405.26_GRCh38_genomic.gtf \
  -o ${dir}/4counts/counts.txt \
  ${dir}/3align/*.bam \
  > ${dir}/4counts/featureCounts.log 2>&1
echo "done" && python3 ${dir}/code/a.py
