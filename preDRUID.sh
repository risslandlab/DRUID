#!/bin/bash
#
# Quantify reads mappings to different genomic features (exons, introns) using
#   the genomecov and intersect tools of bedtools and in-house python script
#   to summarise results.
# This script assumes the existence of three variables prior to running:
#    GTF_NO_EXON - contains the path to a gtf with exon coordinates
#    GTF_INTRON  - contains the path to a gtf with intron coordinates
#    BAM_PATH    - contains the path to a folder with bam files
# This script can be easily modified to take advantage of tools such as gnu
#   parallel to speed up processing.

# set directory for temp file storage, use RAM if possible cause why not
TEMP=/dev/shm/;

# create temp strand specific gtf files
GTF_NO_EXON_NEG=${TEMP}/exon.neg.gtf;
GTF_NO_EXON_POS=${TEMP}/exon.pos.gtf;
cat $GTF_NO_EXON | grep -P '\s-\s' | sort -k1,1 -k4,4n > $GTF_NO_EXON_NEG;
cat $GTF_NO_EXON | grep -P '\s\+\s' | sort -k1,1 -k4,4n > $GTF_NO_EXON_POS;

GTF_INTRON_NEG=${TEMP}/intron.neg.gtf;
GTF_INTRON_POS=${TEMP}/intron.pos.gtf;
cat $GTF_INTRON | grep -P '\s-\s' | sort -k1,1 -k4,4n > $GTF_INTRON_NEG;
cat $GTF_INTRON | grep -P '\s\+\s' | sort -k1,1 -k4,4n > $GTF_INTRON_POS;

mkdir intersect;

# for each BAM file that we have
for BAM in ${BAM_PATH}/*.bam;
do 
	BASE=$(basename $BAM);

	# create temporary bam file containing only unique mappers
	samtools view -h $BAM | \
		grep -P "^@|NH:i:1" | \
		samtools view -bS - > ${TEMP}/${BASE};

	# create a bedgraph file for both positive and negative strands
	bedtools genomecov -ibam ${TEMP}/${BASE} -bg -split -strand - | \
		sort -k1,1 -k2,2n > ${TEMP}/${BASE%bam}neg.bedgraph;
	bedtools genomecov -ibam ${TEMP}/${BASE} -bg -split -strand + | \
		sort -k1,1 -k2,2n > ${TEMP}/${BASE%bam}pos.bedgraph;
	rm ${TEMP}/${BASE}; 

	# determine intersection between bedgraph and gtfs
	bedtools intersect -a $GTF_NO_EXON_NEG \
		-b ${TEMP}/${BASE%bam}pos.bedgraph -wo -sorted \
		-g $GENOME_SIZES > ${TEMP}/${BASE}_exon_pos.intersect
	../../Scripts/coverage.py ${TEMP}/${BASE}_exon_pos.intersect > \
		intersect/${BASE%.bam}_exon_pos.csv;
	rm ${TEMP}/${BASE}_exon_pos.intersect;
	bedtools intersect -a $GTF_NO_EXON_POS \
		-b ${TEMP}/${BASE%bam}neg.bedgraph -wo -sorted \
		-g $GENOME_SIZES > ${TEMP}/${BASE}_exon_neg.intersect
	../../Scripts/coverage.py ${TEMP}/${BASE}_exon_neg.intersect > \
	intersect/${BASE%.bam}_exon_neg.csv;
	rm ${TEMP}/${BASE}_exon_neg.intersect;

	bedtools intersect -a $GTF_INTRON_NEG \
		-b ${TEMP}/${BASE%bam}pos.bedgraph -wo -sorted \
		-g $GENOME_SIZES > ${TEMP}/${BASE}_intron_pos.intersect
	../../Scripts/coverage.py ${TEMP}/${BASE}_intron_pos.intersect > \
		intersect/${BASE%.bam}_intron_pos.csv;
	rm ${TEMP}/${BASE}_intron_pos.intersect;
	bedtools intersect -a $GTF_INTRON_POS \
		-b ${TEMP}/${BASE%bam}neg.bedgraph -wo -sorted \
		-g $GENOME_SIZES > ${TEMP}/${BASE}_intron_neg.intersect
	../../Scripts/coverage.py ${TEMP}/${BASE}_intron_neg.intersect > \
		intersect/${BASE%.bam}_intron_neg.csv;
	rm ${TEMP}/${BASE}_intron_neg.intersect;

	rm ${TEMP}/${BASE%bam}pos.bedgraph ${TEMP}/${BASE%bam}neg.bedgraph;
done

# remove strand specific gtf files
rm $GTF_NO_EXON_NEG $GTF_NO_EXON_POS $GTF_INTRON_NEG $GTF_INTRON_POS
