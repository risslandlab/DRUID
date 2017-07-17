#!/bin/bash

if [ $# -eq 0 ] || [ ! -e $1 ]
then
  echo "Please pass file with path to raw data files as first arg!"
  exit 1 
fi

mkdir TrimmedReads
cp ~/Trimmomatic-0.36/adapters/TruSeq3-SE.fa ./TrimmedReads
while read FILE; do
  # run Trimmomatic on each FILE
  BASE=`basename ${FILE}`;
  echo -ne "Running Trimmomatic on: ${FILE}: ";
  java -jar ~/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 $FILE ./TrimmedReads/$BASE ILLUMINACLIP:./TrimmedReads/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done <$1
