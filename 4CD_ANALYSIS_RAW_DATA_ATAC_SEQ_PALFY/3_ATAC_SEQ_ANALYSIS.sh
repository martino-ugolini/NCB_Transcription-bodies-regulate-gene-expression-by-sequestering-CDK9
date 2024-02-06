#!/bin/bash
#SBATCH --job-name ATAC_Seq
#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 24
#SBATCH --output ATAC_Seq.%A_%a.log
#SBATCH --array 1-10

INDIR_RAW="/scratch/mugolini/ATAC/RAW_DATA/"
OUTDIR_TRIMMOMATIC="/scratch/mugolini/ATAC/Trimmomatic_Output/"
ADAPTER_TRIMMOMATIC="/scratch/mugolini/ATAC/Trimmomatic_adapters/"

OUTDIR_FASTQC="/scratch/mugolini/ATAC/FastQC_Output/"

INDIR_INDEX="/scratch/mugolini/GENOME_ZEBRAFISH/genome_zebrafish"
OUTDIR_BOWTIE2="/scratch/mugolini/ATAC/Bowtie2_Output/"

GENOME_FAI="/scratch/mugolini/GRCz11_105_GENOME/Danio_rerio.GRCz11.105.dna.primary_assembly.fa.fai"

FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${INDIR_RAW}filelist.txt)

module load gcc
module load trimmomatic
module load fastqc
module load bowtie2
module load samtools
module load picard
module load bedtools2

echo "ATAC-Seq Analysis"
echo "Prepare a filelist.txt file in INDIR_RAW with all the file names (without the _1.fastq.gz, _2.fastq.gz)"
echo "Running analysis on $FILE..."

echo "Trimmomatic trimming..."
mkdir -p $OUTDIR_TRIMMOMATIC
mkdir -p ${OUTDIR_TRIMMOMATIC}Singletons/
trimmomatic PE -threads ${SLURM_CPUS_PER_TASK} \
"${INDIR_RAW}${FILE}_1.fastq.gz" "${INDIR_RAW}${FILE}_2.fastq.gz" \
"${OUTDIR_TRIMMOMATIC}${FILE}_1.fastq.gz" "${OUTDIR_TRIMMOMATIC}Singletons/${FILE}_1.singletons.fastq.gz" \
"${OUTDIR_TRIMMOMATIC}${FILE}_2.fastq.gz" "${OUTDIR_TRIMMOMATIC}Singletons/${FILE}_2.singletons.fastq.gz" \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ILLUMINACLIP:"${ADAPTERS_PATH}TruSeq3-PE.fa":2:30:10:2:True

echo "Running FastQC..."
mkdir -p ${OUTDIR_FASTQC}
fastqc -t ${SLURM_CPUS_PER_TASK} -o $OUTDIR_FASTQC ${OUTDIR_TRIMMOMATIC}${FILE}_1.fastq.gz
fastqc -t ${SLURM_CPUS_PER_TASK} -o $OUTDIR_FASTQC ${OUTDIR_TRIMMOMATIC}${FILE}_2.fastq.gz

echo "Mapping paired end reads..."
# Note: --dovetail option is important to not discard PE reads where the mates extend beyond each other
# (which is likely to happen in case of short fragments with long read lengths)
#bowtie2-build  <genome.fa>  <genomeIndexName>
mkdir -p $OUTDIR_BOWTIE2
(bowtie2 --threads 22 -X 2000 --no-mixed --no-discordant --no-unal --dovetail -x $INDIR_INDEX -1 ${OUTDIR_TRIMMOMATIC}${FILE}_1.fastq.gz -2 ${OUTDIR_TRIMMOMATIC}${FILE}_2.fastq.gz) 2> ${OUTDIR_BOWTIE2}${FILE}.alignment.log \
| samtools view -b -h - > ${OUTDIR_BOWTIE2}${FILE}_aligned.bam

echo "Filtering alignments by quality and mitochondrial origin..."
# deduplication and q30 filtering of reads (because everyone does it), remove reads that are PCR or optical duplicates, keep reads that are in proper pair and have high mapping quality >=30
# and filter out mitchondrial alignments -> these affect the global background level for peak calling
samtools view --threads $SLURM_CPUS_PER_TASK -h -q 30 -f 2 ${OUTDIR_BOWTIE2}${FILE}_aligned.bam | egrep -v "MT" | samtools view -h -b -t $GENOME_FAI - > ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.bam
samtools sort -@ $SLURM_CPUS_PER_TASK -o ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.sorted.bam ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.bam
rm ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.bam

echo "Mark and remove duplicated alignments..."
picard MarkDuplicates INPUT=${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.sorted.bam OUTPUT=${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupmarked.bam VALIDATION_STRINGENCY=SILENT METRICS_FILE=metrics.${FILE}.txt 2> ${FILE}.mark.dup

echo "Quality check: The two numbers should be equal"
samtools view --threads $SLURM_CPUS_PER_TASK -c ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.sorted.bam
samtools view --threads $SLURM_CPUS_PER_TASK -c ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupmarked.bam

samtools view --threads $SLURM_CPUS_PER_TASK -F 1024 -h -o ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupremoved.bam ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupmarked.bam
rm ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupmarked.bam

echo "Sorting by name..."
samtools sort -@ $SLURM_CPUS_PER_TASK -n -o ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupremoved.sortedbyname.bam ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupremoved.bam
rm ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupremoved.bam

echo "Calculating statistics..."
rm ${FILE}_statistics.csv
# c1 and c2 are the number of fragments (each fragment is composed of a fwd read in file _1 and a rev read in file_2)
c1=$(zcat ${INDIR_RAW}${FILE}_1.fastq.gz | echo $((`wc -l`/4)))
c2=$(zcat ${INDIR_RAW}${FILE}_2.fastq.gz | echo $((`wc -l`/4)))
echo $c1;echo $c2
if (($c1 == $c2)); then echo "Fastq_raw;$c1" >> "${FILE}_statistics.csv"; else echo "Fastq_raw;ERROR!" >> "${FILE}_statistics.csv"; fi

c1=$(zcat ${OUTDIR_TRIMMOMATIC}${FILE}_1.fastq.gz | echo $((`wc -l`/4)))
c2=$(zcat ${OUTDIR_TRIMMOMATIC}${FILE}_2.fastq.gz | echo $((`wc -l`/4)))
echo $c1;echo $c2
if (($c1 == $c2)); then echo "Fastq_trimmomatic;$c1" >> "${FILE}_statistics.csv"; else echo "Fastq_trimmomatic;ERROR!" >> "${FILE}_statistics.csv"; fi

# c3 and c4 are the number of fragments (each fragment is composed of a fwd read in file _1 and a rev read in file_2)
c3=$(samtools view --threads $SLURM_CPUS_PER_TASK -c ${OUTDIR_BOWTIE2}${FILE}_aligned.bam)
c4=$(samtools view --threads $SLURM_CPUS_PER_TASK -c -F 1024 ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.sorted.bam)
c5=$(samtools view --threads $SLURM_CPUS_PER_TASK -c ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupremoved.sortedbyname.bam)
echo $c3; echo $c4; echo $c5
echo "mapped;$(expr "$c3" / 2)" >> "${FILE}_statistics.csv"
echo "HighQuality & NoMT;$(expr "$c4" / 2)" >> "${FILE}_statistics.csv"
echo "Not duplicated;$(expr "$c5" / 2)" >> "${FILE}_statistics.csv"

echo "Script finished"

# https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/
