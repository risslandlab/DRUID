#!/bin/bash

mkdir BAM;

for BAM in `find STAR -iname "*.bam"`
do
  NAME=`echo $BAM | cut -d/ -f2`
  cp $BAM BAM/${NAME}.bam
done
