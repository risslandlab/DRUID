mkdir htseq_exons_introns;
parallel --no-notice \
'htseq-count \
--format=bam \
--order=pos \
--stranded=reverse \
--minaqual=10 \
--type=exon \
--idattr=gene_id \
--mode=intersection-strict \
--quiet \
BAM/{} \
$GTF_EXON > \
./htseq_exons_introns/{/.}_exon.tsv; \
htseq-count \
--format=bam \
--order=pos \
--stranded=reverse \
--minaqual=10 \
--type=exon \
--idattr=gene_id \
--mode=union \
--quiet \
BAM/{} \
$GTF_INTRON > \
./htseq_exons_introns/{/.}_intron.tsv;' \
::: `ls BAM | grep ".bam$"`
