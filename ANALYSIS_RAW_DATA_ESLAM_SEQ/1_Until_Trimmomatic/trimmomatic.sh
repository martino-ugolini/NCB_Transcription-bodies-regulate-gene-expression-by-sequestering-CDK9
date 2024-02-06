#!/bin/bash

# NOTE: The adapter sequences need to be downloaded from the GitHub website of Trimmomatic and kept in an "Trimmomatic_adapters" folder.

INDIR=$1
OUTDIR=$2
T=$3
FILE=$4
ADAPTERS_PATH=$5

if [[ -d "$INDIR" ]]
then
    echo "Input directory found  on your filesystem."
else
		echo "Input directory could not be found  on your filesystem. The script will be aborted"
		exit 1
fi

mkdir -p $OUTDIR
mkdir -p ${OUTDIR}Singletons/

echo "This script will analyze the raw fastq files in $INDIR using Trimmomatic and perform both adapter- and quality-based trimming"
echo "The output files will be saved in $OUTDIR"

trimmomatic PE -threads $T \
"${INDIR}${FILE}_R1.fastq.gz" "${INDIR}${FILE}_R2.fastq.gz" \
"${OUTDIR}${FILE}_R1.fastq.gz" "${OUTDIR}Singletons/${FILE}_R1.singletons.fastq.gz" \
"${OUTDIR}${FILE}_R2.fastq.gz" "${OUTDIR}Singletons/${FILE}_R2.singletons.fastq.gz" \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ILLUMINACLIP:"${ADAPTERS_PATH}TruSeq3-PE.fa":2:30:10
