#!/bin/bash

INDIR=$1
OUTDIR=$2
T=$3
FILE=$4


if [[ -d "$INDIR" ]]
then
    echo "Input directory found  on your filesystem."
else
		echo "Input directory could not be found  on your filesystem. The script will be aborted"
		exit 1
fi

mkdir -p $OUTDIR

fastqc -t $T -o $OUTDIR ${INDIR}${FILE}_R1.fastq.gz
fastqc -t $T -o $OUTDIR ${INDIR}${FILE}_R2.fastq.gz
