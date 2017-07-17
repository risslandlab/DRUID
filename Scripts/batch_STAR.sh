#!/bin/bash

if [ $# -eq 0 ] || [ ! -e $1 ]
then
  echo "Please pass file with path to raw data files as first arg!"
  exit 1 
fi

mkdir STAR;

while read FILE; do
  # run STAR for each file
  echo -ne "Running STAR on ${FILE}: ";
  BASE=`basename "${FILE}" .fastq.gz`;
  cd STAR;
  mkdir "${BASE}" && cd $_;
    
  STAR \
    --runThreadN 8 \
    --genomeDir ${GENOME} \
    --readFilesIn ${FILE} \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNoverLmax 0.05 \
    --outFilterScoreMinOverLread 0.75 \
    --outFilterMatchNminOverLread 0.85 \
    --alignIntronMax 1 \
    --outFilterIntronMotifs RemoveNoncanonical \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --limitBAMsortRAM 30000000000 \
    --genomeLoad LoadAndKeep;
  cd ../..;
done <$1
